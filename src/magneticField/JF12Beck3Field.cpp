#include "crpropa/magneticField/JF12Beck3Field.h"

#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#include <iostream>

namespace crpropa {


#ifdef CRPROPA_HAVE_FFTW3F

// turbulent field with Kolmogorov spectrum, B_rms = 1 and Lc = 220 parsec, used for iso&aniiso vector fields
void JF12Beck3Field::randomG(int seed, int Lmin, int Lmax) {
	// useStriated = true;
	// useTurbulent = true;
	G = new VectorGrid(Vector3d(0.), 256, 4 * parsec);
	initTurbulence(G, 1, Lmin * parsec, Lmax * parsec, -11./3., seed);
}

#endif


void JF12Beck3Field::setG(ref_ptr<VectorGrid> grid) {
	// useStriated = true;
	// useTurbulent = true;
	G = grid;
}


ref_ptr<VectorGrid> JF12Beck3Field::getG() {
	return G;
}


Vector3d JF12Beck3Field::getStriatedField(const Vector3d& pos) const {
	if ((pos.getR() < 1 * kpc) or (pos.getR() > 20 * kpc))
		return Vector3d(0, 0, 0); // 0 field for d < 1 kpc or d > 20 kpc

  return getRegularField(pos).dot(G->interpolate(pos)) * getRegularField(pos) / getRegularField(pos).getR() * sqrtbeta * sqrt(3./2.);
}


Vector3d JF12Beck3Field::getAniField(const Vector3d& pos) const {
  return getStriatedField(pos);
}


Vector3d JF12Beck3Field::getTurbulentField(const Vector3d& pos) const {
	return (G->interpolate(pos) * getTurbulentStrength(pos));
}


Vector3d JF12Beck3Field::getIsoField(const Vector3d& pos) const {
  return getTurbulentField(pos);
}


Vector3d JF12Beck3Field::getField(const Vector3d& pos) const {
	Vector3d b(0.);
	if (useTurbulent)
		b += 0.6 * getTurbulentField(pos);
	if (useStriated)
		b += 0.3 * getStriatedField(pos);
	else if (useRegular)
		b += getRegularField(pos);
	return b;
}

} // namespace crpropa
