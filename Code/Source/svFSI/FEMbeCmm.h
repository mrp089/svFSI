#pragma once
//=============================================================================
// This plugin example illustrates how to create a new material. 
// It requires FEBio 2.5 (or up)
//
// Author Steve Maas
// Copyright (c) 2015 - 2016
// All rights reserved
//
//=============================================================================

//-----------------------------------------------------------------------------
// We need to include this file since our new material class will inherit from
// FEElasticMaterial which is defined in this include files.
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// This material class implements a neo-Hookean constitutive model. 
// Since it is a (compressible, coupled) hyper-elatic material, it needs to inherit
// from FEElasticMaterial. 
class FEMbeCmm : public FEElasticMaterial
{
public:
	// The constructor is called when an instance of this class is created.
	// All classes registered by the framework must take the FEModel* as the only
	// parameter in the constructor, even if the class does not need it (which most often
	// will be the case). For material classes, the FEModel parameter is passed to the 
	// base class in the initialization list.
	FEMbeCmm(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	// The following three functions are the only functions that have to be defined.
	// (Actually, the Init function is optional). They are defined as virtual member 
	// functions in the base classes and will be called by FEBio to invoke the specific
	// class implementations.

	// This function calculates the spatial (i.e. Cauchy or true) stress.
	// It takes one parameter, the FEMaterialPoint and returns a mat3ds object
	// which is a symmetric second-order tensor.
	virtual mat3ds Stress(FEMaterialPoint& pt);

	// This function calculates the spatial elasticity tangent tensor. 
	// It takes one parameter, the FEMaterialPoint and retursn a tens4ds object
	// which is a fourth-order tensor with major and minor symmetries.
	virtual tens4dss Tangent(FEMaterialPoint& pt);
};
