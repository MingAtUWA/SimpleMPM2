#ifndef __Linear_Elasticity_H__
#define __Linear_Elasticity_H__

#include "ConstitutiveModel.h"

int linear_elasticity_integration_function(ConstitutiveModel *_self, double dstrain[6]);

class LinearElasticity : public ConstitutiveModel
{
friend int linear_elasticity_integration_function(ConstitutiveModel *_self, double dstrain[6]);
protected:
	double E, niu;

public:
	LinearElasticity() :
		ConstitutiveModel(linear_elasticity_integration_function) {}
	~LinearElasticity() {}

	void set_param(double _E, double _niu);
};

#endif