#include "cppm_framework.h"
#include <mitsuba/core/fstream.h>

MTS_NAMESPACE_BEGIN

struct SPPMGatherPoint : public SPPMFrameworkGatherPoint {
    Spectrum flux, sumEmission;
    Float N;
    inline SPPMGatherPoint() : SPPMFrameworkGatherPoint() { }
};

class SPPMFrameworkIntegrator : public SPPMFramework<SPPMGatherPoint> {
public:
    SPPMFrameworkIntegrator(const Properties &props) : SPPMFramework(props) {
        m_initRadiusFile = props.getString("initRadiusFile", "");
        m_hasInitRadiusFile = m_initRadiusFile.length() != 0;
    }
    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        if (m_hasInitRadiusFile) {
            ref<FileStream> filestream = new FileStream(
                (boost::filesystem::current_path() / m_initRadiusFile).c_str(),
                FileStream::EReadOnly);
            ref<Bitmap> bitmap = new Bitmap(Bitmap::EOpenEXR, filestream);
            size_t nEntries = 
                (size_t) bitmap->getSize().x *
                (size_t) bitmap->getSize().y;
            m_initRadiusData = new Float[nEntries];
            switch (bitmap->getComponentFormat()) {
                case Bitmap::EFloat16: { 
                        half *bData = bitmap->getFloat16Data();
                        for (size_t i=0; i<nEntries; ++i) {
                            m_initRadiusData[i] = (Float) (*bData);
                            bData += bitmap->getChannelCount();
                          }
                    }
                    break;
                case Bitmap::EFloat32: {
                        float *bData = bitmap->getFloat32Data();
                        for (size_t i=0; i<nEntries; ++i) {
                            m_initRadiusData[i] = (Float) (*bData);
                            bData += bitmap->getChannelCount();
                          }
                    }
                    break;
                case Bitmap::EUInt32: {
                        uint32_t *bData = bitmap->getUInt32Data();
                        for (size_t i=0; i<nEntries; ++i) {
                            m_initRadiusData[i] = (Float) (*bData);
                            bData += bitmap->getChannelCount();
                        }
                    }
                    break;
                default:
                    Log(EError, "Unsupported component format!");
            }
        }
        return SPPMFramework<SPPMGatherPoint>::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
    }
    void postprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        delete [] m_initRadiusData;
        SPPMFramework<SPPMGatherPoint>::postprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
    }
    SPPMFrameworkIntegrator(Stream *stream, InstanceManager *manager)
     : SPPMFramework(stream, manager) { }
    std::string getName() { return "SPPM"; }
    MTS_DECLARE_CLASS()
private:
    inline void initGatherPoint(SPPMGatherPoint &gp) {
        gp.radius = m_hasInitRadiusFile ? m_initRadiusData[gp.pos.y * m_bitmap->getWidth() + gp.pos.x] : 
                    (m_kNN == 0 ? m_initialRadius : 0.f);
        gp.flux = Spectrum(0.f);
        gp.sumEmission = Spectrum(0.f);
        gp.N = 0;
    }
    inline Spectrum updateGatherPoint(SPPMGatherPoint &gp, Sampler *sampler) {
        Float M, N = gp.N;
        Spectrum flux, contrib;
        gp.sumEmission += gp.emission;
        contrib = gp.sumEmission / m_iteration;

        if (gp.depth != -1) {
            M = (Float) m_photonMap->estimateRadianceRaw(
                gp.its, gp.radius, flux, m_maxDepth == -1 ? INT_MAX : m_maxDepth-gp.depth);
            if (M > 0) {
                Float ratio = (N + m_alpha * M) / (N + M);
                gp.radius = gp.radius * std::sqrt(ratio);

                gp.flux = (gp.flux + gp.weight * flux) * ratio;
                gp.N = N + m_alpha * M;
            }
        }

        if (gp.radius > 0)
            contrib += gp.flux / ((Float) m_totalEmitted * gp.radius * gp.radius * M_PI);
        return contrib;
    }
    bool m_hasInitRadiusFile;
    std::string m_initRadiusFile;
    Float *m_initRadiusData;
};

MTS_IMPLEMENT_CLASS_S(SPPMFrameworkIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(SPPMFrameworkIntegrator, "Stochastic progressive photon mapper");

MTS_NAMESPACE_END
