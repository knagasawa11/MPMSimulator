//
//  constitutive_model.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/28.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

#ifndef constitutive_model_hpp
#define constitutive_model_hpp


#include <tuple>
#include <algorithm>
#include "boost/math/tools/roots.hpp"

#include "mpm_state.h"
#include "../Utils/basisfunctions.h"
#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif
#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif


using Vectord = Vector<double, SET::dim>;
using Vectori = Vector<int, SET::dim>;
using Matrixd = Matrix<double, SET::dim>;
using Matrixi = Matrix<int, SET::dim>;

using MatrixXd = MatrixX<double, SET::dim>;

enum class ConstitutiveModels
{
	Elastic,
	SnowPlastic,
	HerschelBulkley
};


struct ConstitutiveModel {
	
	virtual void set_material_parameter() = 0;
	virtual	void set_material_parameter(const MPM::MaterialParameters& mat_params ,const MPM::SimulationParameters& sim_params) = 0;
	virtual Matrixd cal_particle_stress(const Matrixd& F, const double& J, const Matrixd& b) = 0;
	virtual std::tuple<Matrixd, double, Matrixd> cal_F_placstic_correction(const Matrixd& F, const double& Jp, const Matrixd& b) = 0;
	virtual ~ConstitutiveModel() {}
	
};


struct Elastic : ConstitutiveModel {
	
	double nu;
	double E;
	double hardening;
	
	Elastic ()
	: nu(0.2)
	, E(100.0)
	, hardening(1.0)
	{};
	
	void set_material_parameter() override;
	void set_material_parameter(const MPM::MaterialParameters& mat_params ,const MPM::SimulationParameters& sim_params) override;
	Matrixd cal_particle_stress(const Matrixd& F, const double& Jp, const Matrixd& b) override;
	std::tuple<Matrixd, double, Matrixd> cal_F_placstic_correction(const Matrixd& F, const double& Jp, const Matrixd& b) override;
};

struct SnowPlastic : ConstitutiveModel {
	
	double nu;
	double E;
	double hardening;
	
	SnowPlastic ()
	: nu(0.2)
	, E(100.0)
	, hardening(1.0)
	{};
	
	void set_material_parameter() override;
	void set_material_parameter(const MPM::MaterialParameters& mat_params ,const MPM::SimulationParameters& sim_params) override;
	Matrixd cal_particle_stress(const Matrixd& F, const double& Jp, const Matrixd& b) override;
	std::tuple<Matrixd, double, Matrixd> cal_F_placstic_correction(const Matrixd& F, const double& Jp, const Matrixd& b) override;
	
};

struct HerschelBulkley : ConstitutiveModel {
	
	double mu;
	double kappa;
	double eta;
	double h;
	double sigma_Y;
	double dt;

	
	HerschelBulkley ()
	: mu(100.0)
	, kappa(1000.0)
	, eta(0.8)
	, h(0.8)
	, sigma_Y(0.1)
	{};
	
	void set_material_parameter() override;
	void set_material_parameter(const MPM::MaterialParameters& mat_params ,const MPM::SimulationParameters& sim_params ) override;
	Matrixd cal_particle_stress(const Matrixd& F, const double& Jp, const Matrixd& b) override;
	std::tuple<Matrixd, double, Matrixd> cal_F_placstic_correction(const Matrixd& F, const double& Jp, const Matrixd& b) override;
};


inline void Elastic::set_material_parameter(const MPM::MaterialParameters& mat_params, const MPM::SimulationParameters& sim_params)
{
	std::cout << "Elastic" << std::endl;
	nu 				= mat_params.nu;
	E					= mat_params.E;
	hardening	= mat_params.hardening;
}

inline void Elastic::set_material_parameter()
{
	std::cout << "elastic" << std::endl;
	E 				= SET::MAT_PARAM::E;
	nu				= SET::MAT_PARAM::nu;
	hardening	= SET::MAT_PARAM::hardening;
}


inline Matrixd Elastic::cal_particle_stress(const Matrixd& F, const double& Jp, const Matrixd& b)
{
	//constexpr double nu = 0.2;
	//constexpr double E = 10000.000;
	//constexpr double hardening = 10.0;

	const double mu_0 = E / (2 * (1 + nu));
	const double lambda_0 = E * nu / ((1+nu) * (1 - 2 * nu));
	
	auto e = std::exp(hardening * (1.0 - Jp));
	double mu = mu_0 * e;
	double lambda = lambda_0 * e;
	double J = F.determinant();
	
	const Eigen::JacobiSVD<Matrixd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	const Matrixd R = svd.matrixU() * svd.matrixV().transpose();
	const Matrixd FP = (2 * mu * (F-R) * F.transpose() + lambda * (J-1) * J * Matrixd::Identity()) ;

	return std::move(FP);
}



inline std::tuple<Matrixd, double, Matrixd> Elastic::cal_F_placstic_correction(const Matrixd& F, const double& Jp, const Matrixd& b)
{
	const Eigen::JacobiSVD<Matrixd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	const double pre_J = F.determinant();
	Matrixd sig = svd.singularValues().asDiagonal();

	const auto F_new = svd.matrixU() * sig * svd.matrixV().transpose();
	const double pro_Jp = Jp * pre_J / F_new.determinant();
	const auto b_new = b;
	
	return {F_new, pro_Jp, b_new};
}


inline void SnowPlastic::set_material_parameter(const MPM::MaterialParameters& mat_params, const MPM::SimulationParameters& sim_params)
{
	std::cout << "SnowPlastic" << std::endl;
	nu 				= mat_params.nu;
	E					= mat_params.E;
	hardening	= mat_params.hardening;
}


inline void SnowPlastic::set_material_parameter()
{
	std::cout << "SnowPlastic" << std::endl;
	E 				= SET::MAT_PARAM::E;
	nu				= SET::MAT_PARAM::nu;
	hardening	= SET::MAT_PARAM::hardening;
}

inline Matrixd SnowPlastic::cal_particle_stress(const Matrixd& F, const double& Jp, const Matrixd& b)
{
	//constexpr double nu = 0.2;
	//constexpr double E = 10000.000;
	//constexpr double hardening = 10.0;

	const double mu_0 = E / (2 * (1 + nu));
	const double lambda_0 = E * nu / ((1+nu) * (1 - 2 * nu));
	
	auto e = std::exp(hardening * (1.0 - Jp));

	double mu = mu_0 * e;
	double lambda = lambda_0 * e;
	double J = F.determinant();
	
	const Eigen::JacobiSVD<Matrixd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	
	const Matrixd R = svd.matrixU() * svd.matrixV().transpose();
	const Matrixd FP = (2 * mu * (F-R) * F.transpose() + lambda * (J-1) * J * Matrixd::Identity()) ;

	return std::move(FP);

}


inline std::tuple<Matrixd, double, Matrixd> SnowPlastic::cal_F_placstic_correction(const Matrixd& F, const double& Jp, const Matrixd& b)
{
	const Eigen::JacobiSVD<Matrixd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	const double pre_J = F.determinant();
	Matrixd sig = svd.singularValues().asDiagonal();
	for (int i=0; i<SET::dim; ++i) {
		//svd.singularValues().asDiagonal()[i, i] = 1.0;
		sig(i,i) = std::clamp(sig(i,i), 1.0 - 0.025, 1.0 + 0.0075);

	}
	const auto F_new = svd.matrixU() * sig * svd.matrixV().transpose();
	const double pro_Jp = std::clamp(Jp * pre_J / F_new.determinant(), 0.6, 20.0);
	const auto b_new = b;
	
	return {F_new, pro_Jp, b_new};
}

inline void HerschelBulkley::set_material_parameter(const MPM::MaterialParameters& mat_params, const MPM::SimulationParameters& sim_params)
{
	std::cout << "HerschelBulkley" << std::endl;
	kappa 		= mat_params.kappa;
	mu				= mat_params.mu;
	eta				= mat_params.eta;
	h					= mat_params.h;
	sigma_Y		= mat_params.sigma_Y;
	dt				= sim_params.dt;
}

inline void HerschelBulkley::set_material_parameter()
{
	std::cout << "HerschelBulkley" << std::endl;
	kappa 		= SET::MAT_PARAM::kappa;
	mu				= SET::MAT_PARAM::mu;
	eta				= SET::MAT_PARAM::eta;
	h					= SET::MAT_PARAM::h;
	sigma_Y		= SET::MAT_PARAM::sigma_Y;
}

Matrixd HerschelBulkley::cal_particle_stress(const Matrixd& F, const double& Jp, const Matrixd& b)
{
	
	const double J = F.determinant();
	const auto b_bar = std::pow(J,-2.0/SET::dim)*b;
	const auto dev_b = b_bar - Matrixd::Identity()*b_bar.trace()/SET::dim;
	
	// Compute stress with volume dependent and preserving functions ( Yue (2) )
	const Matrixd FP = 0.5*kappa*(J + 1.0)*(J - 1.0)*Matrixd::Identity() + mu * dev_b ;

	return std::move(FP);

}



//TODO: move Functor into the H-B Constitutive struct.
/*
struct Functor {
	double s_pre;
	Matrixd b_bar;
	
	double mu;
	double eta, h, sigma_Y;
	Functor(double s_pre, Matrixd b_bar,double mu, double eta, double h, double sigma_Y) : s_pre(s_pre), b_bar(b_bar), mu(mu), eta(eta), h(h), sigma_Y(sigma_Y) {}
	double operator()(double s) {
		
		const double h_inv		= 1.0/h;
		//constexpr double sigma		= .1;
		const double sigma_rt = std::sqrt(2.0/SET::dim)*sigma_Y;
		return std::pow(eta,h_inv)*(s - s_pre) + (2.0/SET::dim)*b_bar.trace()*mu*SET::dt*std::pow(s - sigma_rt,h_inv);
		
	}
};
*/
 
std::tuple<Matrixd, double, Matrixd> HerschelBulkley::cal_F_placstic_correction(const Matrixd& F, const double& Jp, const Matrixd& b)
{
	
	const double J = F.determinant();
	const auto b_bar = std::pow(J,-2.0/SET::dim)*b;
	const auto dev_b = b_bar - Matrixd::Identity()*b_bar.trace()/SET::dim;
	const double sigma_rt = std::sqrt(2.0/SET::dim)*sigma_Y;
	const double h_inv		= 1.0/h;
	double s_pro = mu*dev_b.norm();
	assert(J>0);
	
	
	//　If yielding condition is violated, compute plastic flow
	if(mu*dev_b.norm() > sigma_rt){
		
	//define plastic flow function
//	std::function<double(double)> sf = [s_pro, b_bar, &_mu = mu, &_eta = eta, &_h = h, &_sigma_Y = sigma_Y, &_dt = dt](double s){
		std::function<double(double)> sf = [s_pro, b_bar, &_mu = mu, &_eta = eta, &_h_inv = h_inv, &_sigma_rt = sigma_rt, &_dt = dt](double s){
		//const double h_inv		= 1.0/_h;
		//const double sigma_rt = std::sqrt(2.0/SET::dim)*_sigma_Y;
		return std::pow(_eta,_h_inv)*(s - s_pro) + (2.0/SET::dim)*b_bar.trace()*_mu*_dt*std::pow(s - _sigma_rt,_h_inv);
	};

	boost::math::tools::eps_tolerance<double>
			 tol(std::numeric_limits<double>::digits);
	boost::uintmax_t max = 100;
	std::pair<double, double> r = boost::math::tools::bisect(sf,sigma_rt , mu*dev_b.norm(), tol, max);
	s_pro = (r.first + r.second) / 2.0;
	//std::cout << sigma_rt << " " << s_pro << " " << mu*dev_b.norm() << std::endl;
	//std::cout << sf(s_pro) << std::endl;
	}

	
	Matrixd b_bar_new;
	Matrixd b_bar_norm;
	
	//　Compute b_bar from s and normalizing b_bar
	if(s_pro !=0){
	b_bar_new = (1.0/mu) * s_pro * ( dev_b / dev_b.norm() ) + (1.0/SET::dim)*b_bar.trace()*Matrixd::Identity();
	b_bar_norm = std::pow(b_bar_new.determinant(),-1.0/SET::dim)*b_bar_new;
	}else{
	b_bar_new = b_bar;
	b_bar_norm = std::pow(b_bar_new.determinant(),-1.0/SET::dim)*b_bar_new;
	}
	const auto b_new = b_bar_norm*std::pow(J,2.0/SET::dim);
	
	//std::cout << b_new.determinant() << " " << F.determinant()*F.determinant() << " " << std::endl;
	
	return {F, Jp, b_new};
}


#endif /* constitutive_model_h */
