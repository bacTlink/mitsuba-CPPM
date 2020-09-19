/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/cppmgatherproc.h>
#include <mitsuba/render/renderqueue.h>

#include <boost/filesystem.hpp>
#include <boost/thread/xtime.hpp>

#if defined(MTS_OPENMP)
# include <omp.h>
#endif

MTS_NAMESPACE_BEGIN

struct SPPMFrameworkGatherPoint {
    Intersection its;
    Float radius;
    Spectrum weight;
    Spectrum emission;
    int depth;
    Point2i pos;

    inline SPPMFrameworkGatherPoint() : radius(0.f), weight(0.0f), emission(0.0f) { }
};

template <typename GatherPoint>
class SPPMFramework : public Integrator {
public:
    /// Represents one individual PPM gather point including relevant statistics
    SPPMFramework(const Properties &props) : Integrator(props) {
        /* Initial photon query radius (0 = infer based on scene size and sensor resolution) */
        m_initialRadius = props.getFloat("initialRadius", 0);
        m_stepSnapshot = props.getInteger("stepSnapshot", 0);
        m_kNN = props.getInteger("kNN", 0);
        /* Alpha parameter from the paper (influences the speed, at which the photon radius is reduced) */
        m_alpha = props.getFloat("alpha", .66666667);
        /* Number of photons to shoot in each iteration */
        m_photonCount = props.getInteger("photonCount", 250000);
        /* Granularity of the work units used in parallelizing the
           particle tracing task (default: choose automatically). */
        m_granularity = props.getInteger("granularity", 0);
        /* Longest visualized path length (<tt>-1</tt>=infinite). When a positive value is
           specified, it must be greater or equal to <tt>2</tt>, which corresponds to single-bounce
           (direct-only) illumination */
        m_maxDepth = props.getInteger("maxDepth", -1);
        /* Depth to start using russian roulette */
        m_rrDepth = props.getInteger("rrDepth", 3);
        /* Indicates if the gathering steps should be canceled if not enough photons are generated. */
        m_autoCancelGathering = props.getBoolean("autoCancelGathering", true);
        /* Maximum number of passes to render. -1 renders until the process is stopped. */
        m_maxPasses = props.getInteger("maxPasses", -1);
        m_mutex = new Mutex();
        if (m_maxDepth <= 1 && m_maxDepth != -1)
            Log(EError, "Maximum depth must be set to \"2\" or higher!");
        if (m_maxPasses <= 0 && m_maxPasses != -1)
            Log(EError, "Maximum number of Passes must either be set to \"-1\" or \"1\" or higher!");
    }

    SPPMFramework(Stream *stream, InstanceManager *manager)
     : Integrator(stream, manager) { }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Integrator::serialize(stream, manager);
        Log(EError, "Network rendering is not supported!");
    }

    void cancel() {
        m_running = false;
    }

    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
        if (m_initialRadius == 0) {
            /* Guess an initial radius if not provided
              (use scene width / horizontal or vertical pixel count) * 5 */
            Float rad = scene->getBSphere().radius;
            Vector2i filmSize = scene->getSensor()->getFilm()->getSize();

            m_initialRadius = std::min(rad / filmSize.x, rad / filmSize.y) * 100;
        }
        m_tick = boost::posix_time::microsec_clock::local_time();

        boost::filesystem::path path = scene->getDestinationFile();
        std::string prefix = path.filename().string();
        std::string filename = prefix + "_timelog.txt";
        m_timefile = fopen((boost::filesystem::current_path() / filename).string().c_str(), "w");
        return true;
    }

    void postprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        Integrator::postprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
        fclose(m_timefile);
    }

    virtual bool render(Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID, int unused) {
        ref<Scheduler> sched = Scheduler::getInstance();
        ref<Sensor> sensor = scene->getSensor();
        ref<Film> film = sensor->getFilm();
        size_t nCores = sched->getCoreCount();
        Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SSE_STR ") ..",
            film->getCropSize().x, film->getCropSize().y,
            nCores, nCores == 1 ? "core" : "cores");

        Vector2i cropSize = film->getCropSize();

        m_gatherBlocks.clear();
        m_running = true;

        ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), Properties("independent")));

        int blockSize = scene->getBlockSize();

        /* Allocate memory */
        m_bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
        m_bitmap->clear();
        for (int yofs=0; yofs<cropSize.y; yofs += blockSize) {
            for (int xofs=0; xofs<cropSize.x; xofs += blockSize) {
                m_gatherBlocks.push_back(std::vector<GatherPoint>());
                m_offset.push_back(Point2i(xofs, yofs));
                std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[m_gatherBlocks.size()-1];
                int nPixels = std::min(blockSize, cropSize.y-yofs)
                            * std::min(blockSize, cropSize.x-xofs);
                gatherPoints.resize(nPixels);
            }
        }
        initGatherPointPos(scene);

        /* Create a sampler instance for every core */
        std::vector<SerializableObject *> samplers(sched->getCoreCount());
        for (size_t i=0; i<sched->getCoreCount(); ++i) {
            ref<Sampler> clonedSampler = sampler->clone();
            clonedSampler->incRef();
            samplers[i] = clonedSampler.get();
        }

        int samplerResID = sched->registerMultiResource(samplers);

#ifdef MTS_DEBUG_FP
        enableFPExceptions();
#endif

#if defined(MTS_OPENMP)
        Thread::initializeOpenMP(nCores);
#endif

        m_totalEmitted = 0;
        m_totalPhotons = 0;
        m_iteration = 0;
        while (m_running && (m_maxPasses == -1 || m_iteration < m_maxPasses)) {
            ++m_iteration;
            presnapshots(scene, film);
            distributedRTPass(scene, samplers);
            rebuildPhotonMap(job, sceneResID, sensorResID, samplerResID); 
            photonMapPass(samplers, queue, job, film, sceneResID,
                    sensorResID, samplerResID);
            snapshots(scene, film);
        }

#ifdef MTS_DEBUG_FP
        disableFPExceptions();
#endif

        for (size_t i=0; i<samplers.size(); ++i)
            samplers[i]->decRef();

        sched->unregisterResource(samplerResID);
        return true;
    }


    std::string toString() const {
        std::ostringstream oss;
        oss << "SPPMFramework[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  initialRadius = " << m_initialRadius << "," << endl
            << "  alpha = " << m_alpha << "," << endl
            << "  photonCount = " << m_photonCount << "," << endl
            << "  granularity = " << m_granularity << "," << endl
            << "  maxPasses = " << m_maxPasses << endl
            << "]";
        return oss.str();
    }

    virtual std::string getName() = 0;

protected:
    std::string getDestFileName(const Scene *scene, std::string id, int iter) {
        std::stringstream ss;
        boost::filesystem::path path = scene->getDestinationFile();
        std::string prefix = path.filename().string();
        ss << prefix.c_str() << "_" << id << "_"  << iter << "passes";
        std::string filename = ss.str();
        return (boost::filesystem::current_path() / filename).string();
    }

    void snapshot(const Scene *scene, ref<Film> film) {
        film->setDestinationFile(getDestFileName(scene, "render", m_iteration), 0);
        film->develop(scene, 0.f);
    }

    struct GetRadius {
        GetRadius(const Float &x): m_initialRadius(x) {}
        Float operator () (const GatherPoint &gp) {
            return gp.radius;
        }
        Float m_initialRadius;
    };

    template <typename Functor> void saveLogImage(const Scene *scene, const std::string &keyword, Functor functor) {
        ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, m_bitmap->getSize());
        m_bitmap->clear();
        Spectrum *target = (Spectrum *) bitmap->getUInt8Data();
        for (int blockIdx = 0; blockIdx<(int) m_gatherBlocks.size(); ++blockIdx) {
            std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];
            for (size_t i=0; i<gatherPoints.size(); ++i) {
                GatherPoint &gp = gatherPoints[i];
                int p = gp.pos.y * bitmap->getWidth() + gp.pos.x;
                target[p] = Spectrum(functor(gp));
            }
        }
        bitmap->write(getDestFileName(scene, keyword, m_iteration) + ".exr", -1);
    }

    virtual void saveRadius(const Scene *scene) {
        saveLogImage(scene, "bandwidth", GetRadius(m_initialRadius));
    }

    virtual void presnapshots(const Scene *scene, ref<Film> film) {
        boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        fprintf(m_timefile, "%d passes: %ld\n", m_iteration, (now - m_tick).total_milliseconds());
        if (m_iteration != m_maxPasses && m_iteration > 10 &&
                (m_stepSnapshot == 0 || m_iteration % m_stepSnapshot != 0))
            return;
        saveRadius(scene);
    }

    virtual void snapshots(const Scene *scene, ref<Film> film) {
        boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        fprintf(m_timefile, "%d passes: %ld\n", m_iteration, (now - m_tick).total_milliseconds());
        if (m_iteration != m_maxPasses && m_iteration > 10 &&
                (m_stepSnapshot == 0 || m_iteration % m_stepSnapshot != 0))
            return;
        snapshot(scene, film);
    }

    Float determindKNNRadius(Point p) {
        assert(m_kNN > 0);
        typename CPPMPhotonMap::SearchResult result[m_kNN + 1];
        Float sqrRadius = m_initialRadius * m_initialRadius;
        m_photonMap->nnSearch(p, sqrRadius, m_kNN, result);
        return std::sqrt(sqrRadius);
    }

    size_t rebuildPhotonMap(const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        Log(EInfo, "Performing a photon mapping pass %i (" SIZE_T_FMT " photons, " SIZE_T_FMT " particles so far.)",
                m_iteration, m_totalPhotons, m_totalEmitted);
        ref<Scheduler> sched = Scheduler::getInstance();

        /* Generate the global photon map */
        ref<CPPMGatherPhotonProcess> proc = new CPPMGatherPhotonProcess(
            CPPMGatherPhotonProcess::EAllSurfacePhotons, m_photonCount,
            m_granularity, m_maxDepth == -1 ? -1 : m_maxDepth-1, m_rrDepth, true,
            m_autoCancelGathering, job);

        proc->bindResource("scene", sceneResID);
        proc->bindResource("sensor", sensorResID);
        proc->bindResource("sampler", samplerResID);

        sched->schedule(proc);
        sched->wait(proc);

        m_photonMap = proc->getPhotonMap();
        m_photonMap->build();
        m_shotParticles = proc->getShotParticles();
        m_totalEmitted += m_shotParticles;
        Log(EDebug, "Photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: "
            SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());
        return proc->getShotParticles();
    }

    void initGatherPointPos(Scene *scene) {
        ref<Sensor> sensor = scene->getSensor();
        ref<Film> film = sensor->getFilm();
        Vector2i cropSize = film->getCropSize();
        int blockSize = scene->getBlockSize();
        for (int i=0; i<(int) m_gatherBlocks.size(); ++i) {
            std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[i];
            int xofs = m_offset[i].x, yofs = m_offset[i].y;
            int index = 0;
            for (int yofsInt = 0; yofsInt < blockSize; ++yofsInt) {
                if (yofsInt + yofs >= cropSize.y)
                    continue;
                for (int xofsInt = 0; xofsInt < blockSize; ++xofsInt) {
                    if (xofsInt + xofs >= cropSize.x)
                        continue;
                    GatherPoint &gatherPoint = gatherPoints[index++];
                    gatherPoint.pos = Point2i(xofs + xofsInt, yofs + yofsInt);
                    initGatherPoint(gatherPoint);
                }
            }
        }
    }

    void distributedRTPass(Scene *scene, std::vector<SerializableObject *> &samplers) {
        ref<Sensor> sensor = scene->getSensor();
        bool needsApertureSample = sensor->needsApertureSample();
        bool needsTimeSample = sensor->needsTimeSample();
        ref<Film> film = sensor->getFilm();
        Vector2i cropSize = film->getCropSize();
        int blockSize = scene->getBlockSize();

        /* Process the image in parallel using blocks for better memory locality */
        Log(EInfo, "Creating %i gather points", cropSize.x*cropSize.y);
        #if defined(MTS_OPENMP)
            #pragma omp parallel for schedule(dynamic)
        #endif
        for (int i=0; i<(int) m_gatherBlocks.size(); ++i) {
            std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[i];
            #if defined(MTS_OPENMP)
                Sampler *sampler = static_cast<Sampler *>(samplers[mts_omp_get_thread_num()]);
            #else
                Sampler *sampler = static_cast<Sampler *>(samplers[0]);
            #endif

            int xofs = m_offset[i].x, yofs = m_offset[i].y;
            int index = 0;
            for (int yofsInt = 0; yofsInt < blockSize; ++yofsInt) {
                if (yofsInt + yofs >= cropSize.y)
                    continue;
                for (int xofsInt = 0; xofsInt < blockSize; ++xofsInt) {
                    if (xofsInt + xofs >= cropSize.x)
                        continue;
                    Point2 apertureSample, sample;
                    Float timeSample = 0.0f;
                    GatherPoint &gatherPoint = gatherPoints[index++];
                    sampler->generate(gatherPoint.pos);
                    if (needsApertureSample)
                        apertureSample = sampler->next2D();
                    if (needsTimeSample)
                        timeSample = sampler->next1D();
                    sample = sampler->next2D();
                    sample += Vector2((Float) gatherPoint.pos.x, (Float) gatherPoint.pos.y);
                    RayDifferential ray;
                    sensor->sampleRayDifferential(ray, sample, apertureSample, timeSample);
                    Spectrum weight(1.0f);
                    int depth = 1;
                    gatherPoint.emission = Spectrum(0.0f);

                    while (true) {
                        if (scene->rayIntersect(ray, gatherPoint.its)) {
                            if (gatherPoint.its.isEmitter())
                                gatherPoint.emission += weight * gatherPoint.its.Le(-ray.d);

                            if (depth >= m_maxDepth && m_maxDepth != -1) {
                                gatherPoint.depth = -1;
                                break;
                            }

                            const BSDF *bsdf = gatherPoint.its.getBSDF();

                            /* Create hit point if this is a diffuse material or a glossy
                               one, and there has been a previous interaction with
                               a glossy material */
                            if ((bsdf->getType() & BSDF::EAll) == BSDF::EDiffuseReflection ||
                                (bsdf->getType() & BSDF::EAll) == BSDF::EDiffuseTransmission ||
                                (depth + 1 > m_maxDepth && m_maxDepth != -1)) {
                                gatherPoint.weight = weight;
                                gatherPoint.depth = depth;
                                break;
                            } else {
                                /* Recurse for dielectric materials and (specific to SPPM):
                                   recursive "final gathering" for glossy materials */
                                BSDFSamplingRecord bRec(gatherPoint.its, sampler);
                                weight *= bsdf->sample(bRec, sampler->next2D());
                                if (weight.isZero()) {
                                    gatherPoint.depth = -1;
                                    break;
                                }
                                ray = RayDifferential(gatherPoint.its.p,
                                    gatherPoint.its.toWorld(bRec.wo), ray.time);
                                ++depth;
                            }
                        } else {
                            /* Generate an invalid sample */
                            gatherPoint.depth = -1;
                            gatherPoint.emission += weight * scene->evalEnvironment(ray);
                            break;
                        }
                    }
                    sampler->advance();
                }
            }
        }
    }

    virtual void photonMapPass(std::vector<SerializableObject *> &samplers, RenderQueue *queue, const RenderJob *job,
            Film *film, int sceneResID, int sensorResID, int samplerResID) {

        Log(EInfo, "Gathering ..");
        m_totalPhotons += m_photonMap->size();
        film->clear();
        #if defined(MTS_OPENMP)
            #pragma omp parallel for schedule(dynamic)
        #endif
        for (int blockIdx = 0; blockIdx<(int) m_gatherBlocks.size(); ++blockIdx) {
            std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];
            #if defined(MTS_OPENMP)
                Sampler *sampler = static_cast<Sampler *>(samplers[mts_omp_get_thread_num()]);
            #else
                Sampler *sampler = static_cast<Sampler *>(samplers[0]);
            #endif

            Spectrum *target = (Spectrum *) m_bitmap->getUInt8Data();
            for (size_t i=0; i<gatherPoints.size(); ++i) {
                GatherPoint &gp = gatherPoints[i];
                if (gp.depth != -1 && gp.radius == 0.f) gp.radius = determindKNNRadius(gp.its.p);
                target[gp.pos.y * m_bitmap->getWidth() + gp.pos.x] = updateGatherPoint(gp, sampler);
            }
        }
        film->setBitmap(m_bitmap);
        Log(EInfo, "Gathering Done");
        queue->signalRefresh(job);
    }

    virtual void initGatherPoint(GatherPoint &gp) = 0;
    virtual Spectrum updateGatherPoint(GatherPoint &gp, Sampler *sampler) = 0;

    std::vector<std::vector<GatherPoint> > m_gatherBlocks;
    std::vector<Point2i> m_offset;
    ref<Mutex> m_mutex;
    ref<Bitmap> m_bitmap;
    ref<CPPMPhotonMap> m_photonMap;
    Float m_initialRadius, m_alpha;
    int m_photonCount, m_granularity;
    int m_maxDepth, m_rrDepth, m_kNN, m_iteration;
    size_t m_totalEmitted, m_totalPhotons, m_shotParticles;
    bool m_running;
    bool m_autoCancelGathering;
    int m_maxPasses, m_stepSnapshot;
    boost::posix_time::ptime m_tick;
    FILE* m_timefile;
};

MTS_NAMESPACE_END
