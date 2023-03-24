#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Implement MirrorBSDF
        reflect(wo, wi);
        *pdf = 1;
        return reflectance/(abs_cos_theta(*wi));
    }

    void MirrorBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Mirror BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Project 3-2: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        double theta_h = acos(h.z);
        return  exp(-1. * pow(tan(theta_h) / alpha, 2)) / (PI * pow(alpha, 2) * pow(h.z, 4));
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.
        double cos_theta_i = abs_cos_theta(wi);
        // helper variable for better calculations
        Vector3D sum = (eta * eta + k * k);
        Vector3D product = 2. * eta * cos_theta_i;
        double squared = pow(cos_theta_i, 2);
        // literally from website
        Vector3D R_s = (sum - product + squared) / (sum + product + squared);
        Vector3D R_p = (sum * squared - product + 1.) / (sum* squared + product + 1.);
        Vector3D F = (R_s + R_p) / 2.;
        return F;
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Implement microfacet model here.
        // check if wi and wo are valid, if not, return 0
        if (wo.z > 0 && wi.z > 0) {
            Vector3D h = (wo + wi) / (wo + wi).norm();
            Vector3D f = (F(wi) * G(wo, wi) * D(h)) / (4. * wo.z * wi.z);
            return f;
        }
        return Vector3D(0.);

        // return Vector3D();
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        // default cosine sampling code
        // *wi = cosineHemisphereSampler.get_sample(pdf);
        // return MicrofacetBSDF::f(wo, *wi);

        // get sample and compute theta_h/phi_h with r1 and r2
        Vector2D r = sampler.get_sample();
        double r1 = r.x;
        double r2 = r.y;
        double theta_h = atan(sqrt(-1. * pow(alpha, 2) * log(1. - r1)));
        double phi_h = 2. * PI * r2;

        // update wi with the sampled half vector
        Vector3D h(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
        *wi = (2. * dot(h, wo)) * h - wo;
        wi->normalize();

        // check if your sampled wi is valid. If not, return zero pdf and zero BRDF
        if ((*wi).z < 0) {
            *pdf = 0.;
            return Vector3D(0.);
        }
    
        // calculate pdf of sampling h wrt solid angle, define helper variables
        double p_theta_exp = exp(-1. * pow(tan(theta_h) / alpha, 2));
        double p_theta_fraction = (2. * sin(theta_h) / (pow(alpha, 2) * pow(cos(theta_h), 3)));
        double p_theta = p_theta_fraction * p_theta_exp;
        double p_phi = 1. / (2. * PI);
        double p_w_h = p_theta * p_phi / sin(theta_h);

        // calculate final pdf
        *pdf = p_w_h / (4. * dot(*wi, h));

        return MicrofacetBSDF::f(wo, *wi);
    }

    void MicrofacetBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Micofacet BSDF"))
        {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 1
        // Implement RefractionBSDF
        *pdf = 1.;
        bool refracts = refract(wo, wi, ior);
        double eta = (wo.z > 0) ? 1./ior : ior;
        if (refracts) {
            return transmittance / abs_cos_theta(*wi) / (eta * eta);
        }
        return Vector3D();
    }

    void RefractionBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.
        bool refracts = refract(wo, wi, ior);
        if (!refracts) {
            reflect(wo, wi);
            *pdf = 1.;
            return reflectance / abs_cos_theta(*wi);
        }
        else {
            double n1 = 1.;
            double n2 = ior;
            double R_0 = pow((n1 - n2)/(n1 + n2), 2);
            double R_theta = R_0 + (1. - R_0) * pow(1. - abs_cos_theta(wo), 5);
            if (coin_flip(R_theta)) {
                *pdf = R_theta;
                reflect(wo, wi);
                return R_theta * reflectance / abs_cos_theta(*wi);
            }
            *pdf  = 1. - R_theta;
            double eta = (wo.z > 0) ? 1./ior : ior;
            return (1 - R_theta) * transmittance / abs_cos_theta(*wi) / (eta * eta);
        }
        // compute Fresnel coefficient and use it as the probability of reflection
        // - Fundamentals of Computer Graphics page 305
        return Vector3D();
    }

    void GlassBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

        // TODO Project 3-2: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.
        *wi = Vector3D(-wo.x, -wo.y, wo.z);
    }

    bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

        // TODO Project 3-2: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.
        double n = (wo.z > 0) ? 1./ior : ior;
        double oppD = (wo.z > 0) ? -1. : 1.;
        double tir = 1. - (n * n) * (1 - cos_theta(wo) * cos_theta(wo));
        if (tir < 0) {
            return false;
        }
        *wi = Vector3D(-n * wo.x, -n * wo.y, oppD * sqrt(tir));
        return true;

    }

} // namespace CGL
