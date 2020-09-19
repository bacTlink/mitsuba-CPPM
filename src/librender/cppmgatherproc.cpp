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

#include <mitsuba/render/cppmgatherproc.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief This work result implementation stores a sequence of photons, which can be
 * sent over the wire as needed.
 *
 * It is used to implement parallel networked photon tracing passes.
 */
class CPPMPhotonVector : public WorkResult {
public:
    CPPMPhotonVector() { }

    inline void nextParticle() {
        m_particleIndices.push_back((uint32_t) m_photons.size());
    }

    inline void put(const Photon &p) {
        m_photons.push_back(p);
    }

    inline size_t size() const {
        return m_photons.size();
    }

    inline size_t getParticleCount() const {
        return m_particleIndices.size()-1;
    }

    inline size_t getParticleIndex(size_t idx) const {
        return m_particleIndices.at(idx);
    }

    inline void clear() {
        m_photons.clear();
        m_particleIndices.clear();
    }

    inline const Photon &operator[](size_t index) const {
        return m_photons[index];
    }

    void load(Stream *stream) {
        clear();
        size_t count = (size_t) stream->readUInt();
        m_particleIndices.resize(count);
        stream->readUIntArray(&m_particleIndices[0], count);
        count = (size_t) stream->readUInt();
        m_photons.resize(count);
        for (size_t i=0; i<count; ++i)
            m_photons[i] = Photon(stream);
    }

    void save(Stream *stream) const {
        stream->writeUInt((uint32_t) m_particleIndices.size());
        stream->writeUIntArray(&m_particleIndices[0], m_particleIndices.size());
        stream->writeUInt((uint32_t) m_photons.size());
        for (size_t i=0; i<m_photons.size(); ++i)
            m_photons[i].serialize(stream);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "CPPMPhotonVector[size=" << m_photons.size() << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    // Virtual destructor
    virtual ~CPPMPhotonVector() { }
private:
    std::vector<Photon> m_photons;
    std::vector<uint32_t> m_particleIndices;
};

/**
 * This class does the actual photon tracing work
 */
class CPPMGatherPhotonWorker : public ParticleTracer {
public:
    CPPMGatherPhotonWorker(CPPMGatherPhotonProcess::EGatherType type, size_t granularity,
        int maxDepth, int rrDepth) : ParticleTracer(maxDepth, rrDepth, false),
        m_type(type), m_granularity(granularity) { }

    CPPMGatherPhotonWorker(Stream *stream, InstanceManager *manager)
     : ParticleTracer(stream, manager) {
        m_type = (CPPMGatherPhotonProcess::EGatherType) stream->readInt();
        m_granularity = stream->readSize();
    }

    ref<WorkProcessor> clone() const {
        return new CPPMGatherPhotonWorker(m_type, m_granularity, m_maxDepth,
            m_rrDepth);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        ParticleTracer::serialize(stream, manager);
        stream->writeInt(m_type);
        stream->writeSize(m_granularity);
    }

    ref<WorkResult> createWorkResult() const {
        return new CPPMPhotonVector();
    }

    void process(const WorkUnit *workUnit, WorkResult *workResult,
        const bool &stop) {
        m_workResult = static_cast<CPPMPhotonVector *>(workResult);
        m_workResult->clear();
        ParticleTracer::process(workUnit, workResult, stop);
        m_workResult->nextParticle();
        m_workResult = NULL;
    }

    void handleNewParticle() {
        m_workResult->nextParticle();
    }

    void handleSurfaceInteraction(int depth_, int nullInteractions, bool delta,
            const Intersection &its, const Medium *medium,
            const Spectrum &weight) {
        int bsdfType = its.getBSDF()->getType(), depth = depth_ - nullInteractions;
        if (!(bsdfType & BSDF::EDiffuseReflection) && !(bsdfType & BSDF::EGlossyReflection))
            return;

        if ((m_type == CPPMGatherPhotonProcess::ECausticPhotons && depth > 1 && delta)
         || (m_type == CPPMGatherPhotonProcess::ESurfacePhotons && depth > 1 && !delta)
         || (m_type == CPPMGatherPhotonProcess::EAllSurfacePhotons))
            m_workResult->put(Photon(its.p, its.geoFrame.n, -its.toWorld(its.wi), weight, depth));
    }

    void handleMediumInteraction(int depth, int nullInteractions, bool delta,
            const MediumSamplingRecord &mRec, const Medium *medium,
            const Vector &wi, const Spectrum &weight) {
        if (m_type == CPPMGatherPhotonProcess::EVolumePhotons)
            m_workResult->put(Photon(mRec.p, Normal(0.0f, 0.0f, 0.0f),
                -wi, weight, depth-nullInteractions));
    }

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~CPPMGatherPhotonWorker() { }
protected:
    CPPMGatherPhotonProcess::EGatherType m_type;
    size_t m_granularity;
    ref<CPPMPhotonVector> m_workResult;
};

CPPMGatherPhotonProcess::CPPMGatherPhotonProcess(EGatherType type, size_t photonCount,
    size_t granularity, int maxDepth, int rrDepth, bool isLocal, bool autoCancel,
    const void *progressReporterPayload)
    : ParticleProcess(ParticleProcess::ETrace, photonCount, granularity, "Gathering photons",
      progressReporterPayload), m_type(type), m_photonCount(photonCount), m_maxDepth(maxDepth),
      m_rrDepth(rrDepth),  m_isLocal(isLocal), m_autoCancel(autoCancel), m_excess(0), m_numShot(0) {
    m_photonMap = new CPPMPhotonMap(photonCount);
}

bool CPPMGatherPhotonProcess::isLocal() const {
    return m_isLocal;
}

ref<WorkProcessor> CPPMGatherPhotonProcess::createWorkProcessor() const {
    return new CPPMGatherPhotonWorker(m_type, m_granularity, m_maxDepth, m_rrDepth);
}

void CPPMGatherPhotonProcess::processResult(const WorkResult *wr, bool cancelled) {
    if (cancelled)
        return;
    const CPPMPhotonVector &vec = *static_cast<const CPPMPhotonVector *>(wr);
    LockGuard lock(m_resultMutex);

    size_t nParticles = 0;
    for (size_t i=0; i<vec.getParticleCount(); ++i) {
        size_t start = vec.getParticleIndex(i),
               end   = vec.getParticleIndex(i+1);
        ++nParticles;
        for (size_t j=start; j<end; ++j)
            m_photonMap->push_back(vec[j]);
    }
    m_numShot += nParticles;
}

ParallelProcess::EStatus CPPMGatherPhotonProcess::generateWork(WorkUnit *unit, int worker) {
    /* Use the same approach as PBRT for auto canceling */
    LockGuard lock(m_resultMutex);
    if (m_autoCancel && m_numShot > 100000
            && unsuccessful(m_photonCount, m_photonMap->size(), m_numShot)) {
        Log(EInfo, "Not enough photons could be collected, giving up");
        return EFailure;
    }

    return ParticleProcess::generateWork(unit, worker);
}

MTS_IMPLEMENT_CLASS(CPPMGatherPhotonProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS_S(CPPMGatherPhotonWorker, false, ParticleTracer)
MTS_IMPLEMENT_CLASS(CPPMPhotonVector, false, WorkResult)
MTS_NAMESPACE_END
