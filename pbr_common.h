#ifndef __PBR_COMMON_H__
#define __PBR_COMMON_H__

struct pbr_param_t
{
	vector3_t n; // world normal
	vector3_t l; // input light direction
	vector3_t h; // half vector h = normalize((v + l) / 2)
	vector3_t v; // view direction
	float NoL;
	float NoH;
	float NoV;
	float HoV;
	float metallic;
	float roughness;
	vector3_t f0;
};

// Point lobe in off-specular peak direction
vector3_t getoff_specular_peak_reflection_dir(vector3_t normal, vector3_t reflection_vector, float roughness)
{
	float a = roughness * roughness;
	return lerp(normal, reflection_vector, (1 - a) * (sqrt(1 - a) + a));
}

//
vector3_t F_fresenl_schlick(float HoV, const vector3_t& f0)
{
	float e = pow(1.0f - HoV, 5.0f);
	return f0 + (vector3_t::one() - f0) * e;
}

vector3_t F_fresenl_schlick_roughness(float NoV, const vector3_t& f0, float roughness)
{
	float e = pow(1.0f - NoV, 5.0f);
	vector3_t r = vector3_t::one() * (1.0f - roughness);
	r = maxv(r, f0);
	return f0 + (r - f0) * e;
}

// GGX / Trowbridge-Reitz
// [Walter et al. 2007, "Microfacet models for refraction through rough surfaces"]
float D_Trowbridge_Reitz_GGX(float a, float NoH)
{
	float a2 = a * a;
	float d = NoH * NoH * (a2 - 1.0f) + 1.0f;
	return a2 / (cPI * d * d);
}

// Tuned to match behavior of Vis_Smith
// [Schlick 1994, "An Inexpensive BRDF Model for Physically-Based Rendering"]
float V_Schlick_GGX(float roughness, float NoV, float NoL)
{
	// V = G / (NoL * NoV)
	roughness = (roughness + 1.0f) / 2.0f; // for direct shading
	float k = roughness * roughness * 0.5f;
	float Vis_SchlickV = NoV * (1.0f - k) + k;
	float Vis_SchlickL = NoL * (1.0f - k) + k;
	return 0.25f / (Vis_SchlickV * Vis_SchlickL);
}

// Smith term for GGX
// [Smith 1967, "Geometrical shadowing of a random rough surface"]
float V_Smith_GGX(float a, float NoV, float NoL)
{
	// V = G / (NoL * NoV)
	float a2 = a * a;
	float Vis_SmithV = NoV + sqrt(NoV * (NoV - NoV * a2) + a2);
	float Vis_SmithL = NoL + sqrt(NoL * (NoL - NoL * a2) + a2);
	float m = Vis_SmithV * Vis_SmithL;
	return m != 0 ? 1.0f / m : FLT_MAX;
}

// Appoximation of joint Smith term for GGX
// [Heitz 2014, "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"]
float V_SmithJointApprox_GGX(float a, float NoV, float NoL)
{
	float Vis_SmithV = NoL * (NoV * (1 - a) + a);
	float Vis_SmithL = NoV * (NoL * (1 - a) + a);
	return 0.5f / (Vis_SmithV + Vis_SmithL);
}

// [Heitz 2014, "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"]
float V_SmithJoint_GGX(float a, float NoV, float NoL)
{
	float a2 = a * a;
	float Vis_SmithV = NoL * sqrt(NoV * (NoV - NoV * a2) + a2);
	float Vis_SmithL = NoV * sqrt(NoL * (NoL - NoL * a2) + a2);
	return 0.5f / (Vis_SmithV + Vis_SmithL);
}

#endif //__PBR_COMMON_H__

