#include "cppm_framework.h"

// null hypothesis: generalized odd function

MTS_NAMESPACE_BEGIN

namespace {

#define SEC_U 2
#define SEC_V 6

using std::min;
using std::max;

struct CPPMGatherPoint: public SPPMFrameworkGatherPoint {
    size_t photonCount[SEC_U * SEC_V];
    Spectrum sumI, sumEmission;
    Float minPhotonCount;
    size_t totalPhotonCount;
    inline CPPMGatherPoint() : SPPMFrameworkGatherPoint(),
        sumI(0.f), sumEmission(0.f), totalPhotonCount(0) {
        memset(photonCount, 0, sizeof(photonCount));
    }
};

struct CPPMQuery {
    CPPMQuery(CPPMGatherPoint &gp, int maxDepth, Sampler *sampler)
        : photonCount(gp.photonCount), 
        its(gp.its), n(gp.its.shFrame.n), position(gp.its.p),
        weight(gp.weight), squaredRadius(gp.radius * gp.radius), maxDepth(maxDepth), sumFlux(0.f),
        sampler(sampler) {
        bsdf = its.getBSDF();
    }
    inline int getSection(const Point &photonPos) {
        Vector vector = its.toLocal(photonPos - position);
        Float angle = atan2(vector.y, vector.x);
        angle += vector.y < 0 || (vector.y == 0 && vector.x < 0) ? 2 * M_PI : 0;
        int i = max(0, min(SEC_U - 1, (int)(SEC_U * vector.lengthSquared() / squaredRadius)));
        int j = (int)(SEC_V * angle / (M_PI * 2)) % SEC_V;
        return i * SEC_V + j;
    }
    inline int getOppositeSection(const Point &photonPos) {
        Vector vector = its.toLocal(photonPos - position);
        Float angle = atan2(vector.y, vector.x);
        angle += vector.y < 0 ? 2 * M_PI : 0;
        angle += M_PI;
        if (angle >= 2 * M_PI)
            angle -= 2 * M_PI;
        int i = max(0, min(SEC_U - 1, (int)(SEC_U * vector.lengthSquared() / squaredRadius)));
        int j = (int)(SEC_V * angle / (M_PI * 2)) % SEC_V;
        return i * SEC_V + j;
    }
    inline void operator()(const Photon &photon) {
        Normal photonNormal(photon.getNormal());
        Vector wi = -photon.getDirection();
        Float wiDotGeoN = absDot(photonNormal, wi);
        if (photon.getDepth() > maxDepth
            || dot(photonNormal, n) < 1e-1f
            || wiDotGeoN < 1e-2f)
            return;

        BSDFSamplingRecord bRec(its, its.toLocal(wi), its.wi, EImportance);

        Spectrum value = photon.getPower() * bsdf->eval(bRec);
        if (value.isZero())
            return;
        /* Account for non-symmetry due to shading normals */
        value *= std::abs(Frame::cosTheta(bRec.wi) /
            (wiDotGeoN * Frame::cosTheta(bRec.wo)));

        int sec = sampler->next1D() < 0.5 ? getSection(photon.position)
                                          : getOppositeSection(photon.position);
        ++photonCount[sec];

        value *= weight;

        sumFlux += value;
    }
    size_t *photonCount;
    const Intersection &its;
    const Normal &n;
    const Point &position;
    Spectrum weight;
    const BSDF *bsdf;
    Float squaredRadius;
    int maxDepth;
    Spectrum sumFlux;
    Sampler *sampler;
};

const double chi2_90[] = {100, 2.70554, 4.60517, 6.25139, 7.77944, 9.23636, 10.64464, 12.01704, 13.36157, 14.68366, 15.98718, 17.27501, 18.54935, 19.81193, 21.06414, 22.30713, 23.54183, 24.76904, 25.98942, 27.20357, 28.41198, 29.61509, 30.81328, 32.00690, 33.19624, 34.38159, 35.56317, 36.74122, 37.91592, 39.08747, 40.25602, 41.42174, 42.58475, 43.74518, 44.90316, 46.05879, 47.21217, 48.36341, 49.51258, 50.65977, 51.80506, 52.94851, 54.09020, 55.23019, 56.36854, 57.50530, 58.64054, 59.77429, 60.90661, 62.03754};
const double chi2_95[] = {100, 3.84146, 5.99146, 7.81473, 9.48773, 11.07050, 12.59159, 14.06714, 15.50731, 16.91898, 18.30704, 19.67514, 21.02607, 22.36203, 23.68479, 24.99579, 26.29623, 27.58711, 28.86930, 30.14353, 31.41043, 32.67057, 33.92444, 35.17246, 36.41503, 37.65248, 38.88514, 40.11327, 41.33714, 42.55697, 43.77297, 44.98534, 46.19426, 47.39988, 48.60237, 49.80185, 50.99846, 52.19232, 53.38354, 54.57223, 55.75848, 56.94239, 58.12404, 59.30351, 60.48089, 61.65623, 62.82962, 64.00111, 65.17077, 66.33865};
const double chi2_99[] = {100, 6.63490, 9.21034, 11.34487, 13.27670, 15.08627, 16.81189, 18.47531, 20.09024, 21.66599, 23.20925, 24.72497, 26.21697, 27.68825, 29.14124, 30.57791, 31.99993, 33.40866, 34.80531, 36.19087, 37.56623, 38.93217, 40.28936, 41.63840, 42.97982, 44.31410, 45.64168, 46.96294, 48.27824, 49.58788, 50.89218, 52.19139, 53.48577, 54.77554, 56.06091, 57.34207, 58.61921, 59.89250, 61.16209, 62.42812, 63.69074, 64.95007, 66.20624, 67.45935, 68.70951, 69.95683, 71.20140, 72.44331, 73.68264, 74.91947};

} // namespace

class CPPMIntegrator : public SPPMFramework<CPPMGatherPoint> {
public:
    CPPMIntegrator(const Properties &props) : SPPMFramework(props) {
        m_k = props.getFloat("k", 0.8);
        m_beta = props.getFloat("beta", 1.5625);
    }

    CPPMIntegrator(Stream *stream, InstanceManager *manager)
     : SPPMFramework(stream, manager) { }

    std::string getName() { return "CPPM-GOF"; }

    MTS_DECLARE_CLASS()

private:
    inline void initGatherPoint(CPPMGatherPoint &gp) {
        gp.radius = m_kNN == 0 ? m_initialRadius : 0.f;
        gp.minPhotonCount = m_kNN;
    }
    void chiSquaredTest(size_t *photonCount, bool passed[SEC_U]) {
        Float sumO = 0.f, sumSqrO = 0.f, V = 0.f;
        for (int i = 0; i < SEC_U; ++i) {
            for (int j = 0; j < SEC_V; ++j) {
                size_t O = photonCount[i * SEC_V + j];
                sumO += O;
                sumSqrO += O * O;
            }
            int secs = (i + 1) * SEC_V;
            if (sumO == 0)
                passed[i] = true;
            else {
                V = sumSqrO * secs / (double)sumO - sumO;
                passed[i] = V < chi2_95[secs - 1];
            }
        }
    }
    Spectrum updateGatherPoint(CPPMGatherPoint &gp, Sampler *sampler) {
        if (gp.depth != -1) {
            int photonMaxDepth = m_maxDepth == -1 ? INT_MAX : m_maxDepth - gp.depth;
            sampler->generate(gp.pos);
            CPPMQuery query(gp, photonMaxDepth, sampler);
            gp.totalPhotonCount += m_photonMap->executeQuery(gp.its.p, gp.radius, query);
            gp.sumI += query.sumFlux / (M_PI * gp.radius * gp.radius);
            if (gp.totalPhotonCount >= gp.minPhotonCount) {
                bool passed[SEC_U];
                chiSquaredTest(gp.photonCount, passed);
                if (!passed[SEC_U - 1]) {
                    memset(gp.photonCount, 0, sizeof(gp.photonCount));
                    gp.totalPhotonCount = 0;
                    gp.minPhotonCount *= m_beta;
                    bool found = false;
                    for (int i = SEC_U - 1; i >= 0; --i)
                        if (passed[i]) {
                            gp.radius *= std::sqrt((Float)(i + 1) / SEC_U);
                            found = true;
                            break;
                        }
                    if (!found) {
                        gp.radius *= std::sqrt(m_k);
                    }
                }
            }
        }
        gp.sumEmission += gp.emission;
        return gp.sumEmission / m_iteration + gp.sumI / m_totalEmitted;
    }
    Float m_k, m_beta;
};

MTS_IMPLEMENT_CLASS_S(CPPMIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(CPPMIntegrator, "Adapative progressive photon mapper");

MTS_NAMESPACE_END
