#ifndef CRPROPA_JF12BECK3FIELD_H
#define CRPROPA_JF12BECK3FIELD_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/Grid.h"

namespace crpropa {

/**
 @class JF12Beck3
 @brief JF12Beck3 galactic magnetic field model

 Implements the Beck3 magnetic field model, consisting of a large-scale regular (JF12)
 and random (striated) field and a small-scale random (turbulent) field.
 See:
 Beck et al. 2016 New constraints on modelling the random magnetic field of the MW

 The field is defined in the usual galactocentric coordinate system with the
 Galactic center at the origin, the x-axis pointing in the opposite direction of
 the Sun, and the z-axis pointing towards Galactic north.
 */
class JF12Beck3Field: public JF12Field {

protected:
	ref_ptr<VectorGrid> G;

public:

	#ifdef CRPROPA_HAVE_FFTW3F
	  void randomG(int seed, int Lmin = 40, int Lmax = 977.8);
		// Create a random realization for the turbulent field
	#endif

	void setG(ref_ptr<VectorGrid> grid);
	ref_ptr<VectorGrid> getG();

	// Regular and striated field component
	Vector3d getStriatedField(const Vector3d& pos) const;  // both same function
	Vector3d getAniField(const Vector3d& pos) const;

	// Turbulent field component
	Vector3d getTurbulentField(const Vector3d& pos) const; // same functions
	Vector3d getIsoField(const Vector3d& pos) const;

	// All set field components
	Vector3d getField(const Vector3d& pos) const;


};

} // namespace crpropa

#endif // CRPROPA_JF12BECK3FIELD_H
