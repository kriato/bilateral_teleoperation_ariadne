#include <iostream>
#include <fstream>

#include "ariadne.hpp"
#include "yaml-cpp/yaml.h"

#define PI Ariadne::pi

// CHANGE THE DEFINE TO CHANGE THE TELEOPERATION INPUT
#define SINE
// #define STEP 

template<typename T>
void print(T in, bool newline = true)
{
	if (newline)
		std::cout << in << std::endl;
	else
		std::cout << in;
}

#if defined(SINE)
	std::string config_path = "../../config_sine.yml";
	double FREQ, AMP;
#elif defined(STEP)
	std::string config_path = "../../config_step.yml";
	double START_TIME, UP_TIME, AMP;
#endif

// OPERATOR
double PO, DO;
// ENVIRONMENT
double QE, KE, BE;
// SLAVE
double ARM_S, JS, BS, BhS, JhS, KT2CS, KC2VS, PS, DS;
// MASTER
double ARM_M, JM, BM, BhM, JhM, KT2CM, KC2VM, PM, DM;
// SIMULATION
Ariadne::Int TF_CONTINUOUS, TF_DISCRETE;
double TS, STEP_SIZE;
// PLOTS
double Y_MIN, Y_MAX;
// TWO LAYER
double HMAX, HMIN, HD, ALPHA, BETA;
// LOG FILENAME
std::string LOG_FILE;

void load_settings()
{
	YAML::Node config = YAML::LoadFile(config_path);

	QE = config["QE"].as<double>(); 
	KE = config["KE"].as<double>(); 
	BE = config["BE"].as<double>(); 
	PO = config["PO"].as<double>();
	DO = config["DO"].as<double>();

	#if defined(SINE) 
	FREQ = config["FREQ"].as<double>(); 
	AMP = config["AMP"].as<double>();
	#elif defined(STEP) 
	UP_TIME = config["UP_TIME"].as<double>(); 
	START_TIME = config["START_TIME"].as<double>(); 
	AMP = config["AMP"].as<double>();
	#endif	

	ARM_S = config["ARM_S"].as<double>();
	JS = config["JS"].as<double>(); 
	BS = config["BS"].as<double>(); 
	BhS = config["BhS"].as<double>(); 
	JhS = config["JhS"].as<double>();
	KT2CS = config["KT2CS"].as<double>(); 
	KC2VS = config["KC2VS"].as<double>();
	PS = config["PS"].as<double>(); 
	DS = config["DS"].as<double>(); 
	ARM_M = config["ARM_M"].as<double>();
	JM = config["JM"].as<double>(); 
	BM = config["BM"].as<double>(); 
	BhM = config["BhM"].as<double>(); 
	JhM = config["JhM"].as<double>();
	KT2CM = config["KT2CM"].as<double>(); 
	KC2VM = config["KC2VM"].as<double>();
	PM = config["PM"].as<double>(); 
	DM = config["DM"].as<double>();
	TF_CONTINUOUS = config["TF_CONTINUOUS"].as<Ariadne::Int>(); 
	TF_DISCRETE = config["TF_DISCRETE"].as<Ariadne::Int>(); 
	Y_MIN = config["Y_MIN"].as<double>(); 
	Y_MAX = config["Y_MAX"].as<double>();
	TS = config["TS"].as<double>(); 
	STEP_SIZE = config["STEP_SIZE"].as<double>();

	HMAX = config["HMAX"].as<double>();
	HMIN = config["HMIN"].as<double>();
	HD = config["HD"].as<double>();
	ALPHA = config["ALPHA"].as<double>();
	BETA = config["BETA"].as<double>();

	LOG_FILE = config["LOG_FILE"].as<std::string>();

	if (config["DEBUG"].as<bool>()) {
		print("ENVIRONMENT");
		print("QE: " + std::to_string(QE));
		print("KE: " + std::to_string(KE));
		print("BE: " + std::to_string(BE));
		print("\nOPERATOR");
		print("PO: " + std::to_string(PO));
		print("DO: " + std::to_string(DO));

		#if defined(SINE) 
		print("FREQ: " + std::to_string(FREQ));
		print("AMP: " + std::to_string(AMP));
		#elif defined(STEP)
		print("UP_TIME: " + std::to_string(UP_TIME));
		print("START_TIME: " + std::to_string(START_TIME));
		print("AMP: " + std::to_string(AMP)); 
		#endif

		print("\nSLAVE");
		print("ARM_S: " + std::to_string(ARM_S));
		print("JS: " + std::to_string(JS));
		print("BS: " + std::to_string(BS));
		print("BhS: " + std::to_string(BhS));
		print("JhS: " + std::to_string(JhS));
		print("KT2CS: " + std::to_string(KT2CS));
		print("KC2VS: " + std::to_string(KC2VS));
		print("PS: " + std::to_string(PS));
		print("DS: " + std::to_string(DS));
		print("\nMASTER");
		print("ARM_M: " + std::to_string(ARM_M));
		print("JM: " + std::to_string(JM));
		print("BM: " + std::to_string(BM));
		print("BhM: " + std::to_string(BhM));
		print("JhM: " + std::to_string(JhM));
		print("KT2CM: " + std::to_string(KT2CM));
		print("KC2VM: " + std::to_string(KC2VM));
		print("PM: " + std::to_string(PM));
		print("DM: " + std::to_string(DM));
		print("\nSIMULATION");
		print("TS: " + std::to_string(TS));
		print("STEP_SIZE: " + std::to_string(STEP_SIZE));
		print("TF_CONTINUOUS: " + std::to_string(TF_CONTINUOUS));
		print("TF_DISCRETE: " + std::to_string(TF_DISCRETE));
		print("\nPLOTS");
		print("Y_MIN: " + std::to_string(Y_MIN));
		print("Y_MAX: " + std::to_string(Y_MAX));
		print("\nTWO LAYER");
		print("HMAX: " + std::to_string(HMAX));
		print("HMIN: " + std::to_string(HMIN));
		print("HD: " + std::to_string(HD));
		print("ALPHA: " + std::to_string(ALPHA));
		print("BETA: " + std::to_string(BETA));
		print("");
	}
}

// Mimic the overriden operator outside the namespace 
Ariadne::Pair<Ariadne::StringVariable, Ariadne::String> getPair(Ariadne::StringVariable base, std::string location) {{ using namespace Ariadne; return base|location;};}
Ariadne::Pair<Ariadne::StringVariable, Ariadne::String> getPair(std::string base, std::string location) {{ using namespace Ariadne; return Ariadne::StringVariable(base)|Ariadne::StringConstant(location);};}

Ariadne::HybridAutomaton Clock()
{
	Ariadne::RealConstant ts("ts", Ariadne::Decimal(TS));
	Ariadne::RealVariable counter("cnt");

	Ariadne::HybridAutomaton clock("clock");
	Ariadne::DiscreteLocation loc;

	Ariadne::DiscreteEvent clock_event("clock_event");

	clock.new_mode(loc, Ariadne::dot({counter}) = {1});

	clock.new_transition(clock_event, {Ariadne::next(counter)=0}, counter >= ts, Ariadne::EventKind::PERMISSIVE);

	return clock;
}

Ariadne::HybridAutomaton Holder()
{	
	Ariadne::RealVariable position_master_d("qm_d");
	Ariadne::RealVariable position_master("qm");
	Ariadne::RealVariable position_slave_d("qs_d");
	Ariadne::RealVariable position_slave("qs");
	Ariadne::RealVariable velocity_master_d("qm_dot_d");
	Ariadne::RealVariable velocity_master("qm_dot");
	Ariadne::RealVariable velocity_slave_d("qs_dot_d");
	Ariadne::RealVariable velocity_slave("qs_dot");
	Ariadne::RealVariable env_force("fe");
	Ariadne::RealVariable env_force_d("fe_d");

	Ariadne::HybridAutomaton holder("holder");
	Ariadne::DiscreteLocation loc;

	Ariadne::DiscreteEvent clock_event("clock_event");

	holder.new_mode(loc, Ariadne::dot({position_master_d, position_slave_d, velocity_master_d, velocity_slave_d, env_force_d}) = {0, 0, 0, 0, 0});

	holder.new_transition(loc, clock_event, loc, Ariadne::next(
		{position_master_d, position_slave_d, velocity_master_d, velocity_slave_d, env_force_d}) =
		{position_master, position_slave, velocity_master, velocity_slave, env_force});

	return holder;
}

Ariadne::HybridAutomaton OneDelay()
{	
	Ariadne::RealVariable position_master_d("qm_d");
	Ariadne::RealVariable position_slave_d("qs_d");	
	Ariadne::RealVariable position_master_d_prev("qm_d_prev");
	Ariadne::RealVariable position_slave_d_prev("qs_d_prev");

	Ariadne::HybridAutomaton one_delay("one_delay");
	Ariadne::DiscreteLocation loc;

	Ariadne::DiscreteEvent clock_event("clock_event");

	one_delay.new_mode(loc, Ariadne::dot({position_slave_d_prev, position_master_d_prev}) = {0, 0});
	
	one_delay.new_transition(loc, clock_event, loc, Ariadne::next({position_master_d_prev, position_slave_d_prev})={position_master_d, position_slave_d});

	return one_delay;
}

Ariadne::HybridAutomaton PLM()
{	
	Ariadne::RealVariable position_master_d("qm_d");
	Ariadne::RealVariable velocity_master_d("qm_dot");
	Ariadne::RealVariable position_master_d_prev("qm_d_prev");
	Ariadne::RealVariable tau_tl("tau_tlm");
	Ariadne::RealVariable tau_pl("tau_plm");
	Ariadne::RealVariable tau_tlc("tau_tlc");
	Ariadne::RealVariable deltaH_debug("deltaH_m");

	// Ariadne::RealConstant ts("ts", Ariadne::Decimal(TS));
	// Ariadne::RealVariable counter("cnt");

	Ariadne::RealConstant zero("zero", Ariadne::Decimal(0));
	Ariadne::RealConstant minus_one("minus_one", Ariadne::Decimal(-1.0));
	
	Ariadne::RealConstant Hmax("Hmax_m", Ariadne::Decimal(HMAX));
	Ariadne::RealConstant Hmin("Hmin_m", Ariadne::Decimal(HMIN));
	Ariadne::RealConstant Hd("Hd_m", Ariadne::Decimal(HD));
	Ariadne::RealConstant beta("beta_m", Ariadne::Decimal(BETA));
	Ariadne::RealConstant alpha("alpha_m", Ariadne::Decimal(ALPHA));

	Ariadne::RealVariable H("Hm");
	Ariadne::RealVariable H_in("H+m");
	Ariadne::RealVariable H_out("H-m");

	Ariadne::StringVariable plm_name("plm");
	Ariadne::HybridAutomaton plm(plm_name.name());

	int n_locations = 12;
	int n_events = 11;

	Ariadne::DiscreteLocation loc[n_locations];
	for (size_t i = 0; i < n_locations; i++)
	{
		loc[i] = Ariadne::DiscreteLocation(getPair(plm_name,"s" + std::to_string(i)));
	}
	
	Ariadne::DiscreteEvent clock_event("clock_event");
	Ariadne::DiscreteEvent events[n_events];
	for (size_t i = 0; i < n_events; i++)
	{
		events[i] = Ariadne::DiscreteEvent("plm_t" + std::to_string(i));
	}
	
	Ariadne::RealExpression newH = H + H_in;
	Ariadne::RealExpression deltaQ = (position_master_d - position_master_d_prev);
	Ariadne::RealExpression deltaH = deltaQ * tau_tl;
	Ariadne::RealExpression Htilde = newH - deltaH;
	Ariadne::RealExpression tau_tlc_exp = -alpha * (Hd - H) * velocity_master_d;
	Ariadne::RealExpression mod_tau_tl = newH - Hmin / deltaQ;

	Ariadne::RealExpression Htrue = Htilde - beta * Htilde;
	Ariadne::RealExpression Hfalse = zero;
	 
	for (size_t i = 0; i < n_locations; i++)
	{
		plm.new_mode(loc[i], Ariadne::dot({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {0, 0, 0, 0, 0});
	} 

	plm.new_transition(loc[0], events[0], loc[1], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (Htilde > Hmin), Ariadne::EventKind::PERMISSIVE);
	plm.new_transition(loc[0], events[1], loc[2], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (Htilde <= Hmin), Ariadne::EventKind::PERMISSIVE);

	plm.new_transition(loc[1], events[2], loc[3], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (Htrue < Hd), Ariadne::EventKind::PERMISSIVE);
	plm.new_transition(loc[1], events[3], loc[4], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (Htrue >= Hd), Ariadne::EventKind::PERMISSIVE);

	plm.new_transition(loc[4], events[4], loc[5], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE);
	plm.new_transition(loc[4], events[5], loc[6], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (Htrue <= Hmax), Ariadne::EventKind::PERMISSIVE);

	plm.new_transition(loc[2], events[6], loc[7], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (deltaQ == zero), Ariadne::EventKind::PERMISSIVE);
	plm.new_transition(loc[2], events[7], loc[8], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (deltaQ != zero), Ariadne::EventKind::PERMISSIVE);

	plm.new_transition(loc[8], events[8], loc[9], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (tau_tl > zero), Ariadne::EventKind::PERMISSIVE);
	plm.new_transition(loc[8], events[9], loc[10], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (Ariadne::RealExpression(tau_tl) == zero), Ariadne::EventKind::PERMISSIVE);
	plm.new_transition(loc[8], events[10], loc[11], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug}) = {tau_pl, H, H_out, tau_tlc, deltaH_debug}, (tau_tl < zero), Ariadne::EventKind::PERMISSIVE);

	// Clock transitions
	plm.new_transition(loc[3], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug})={tau_tl + tau_tlc_exp, Htrue, beta*Htilde, tau_tlc_exp, deltaH}); // Harvesting ON (no saturation possible)
	plm.new_transition(loc[5], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug})={tau_tl + zero, Hmax, beta*Htilde, zero, deltaH}); // Saturation ON
	plm.new_transition(loc[6], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug})={tau_tl + zero, Htrue, beta*Htilde, zero, deltaH}); // Saturation OFF
	plm.new_transition(loc[7], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug})={zero, Hmin, zero, zero, deltaH}); // No energy available and no movement 
	plm.new_transition(loc[9], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug})={abs(mod_tau_tl), Hmin, zero, zero, deltaH}); // Reducing the torque, sign was positive
	plm.new_transition(loc[10], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug})={zero, Hmin, zero, zero, deltaH}); // tau_tl was zero, therefore sign is zero
	plm.new_transition(loc[11], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, tau_tlc, deltaH_debug})={minus_one * abs(mod_tau_tl), Hmin, zero, zero, deltaH}); // Reducing the torque, sign was negative


	// plm.new_transition(events[0], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde > Hmin) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // T_TT  NOT GONNA HAPPEN!
	// plm.new_transition(events[1], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Htrue, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde > Hmin) && (Htrue < Hd) && (Htrue <= Hmax), Ariadne::EventKind::PERMISSIVE); // T_TF
	// plm.new_transition(events[2], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + zero, Hmax, beta*Htilde, zero}, (counter >= ts) && (Htilde > Hmin) && (Htrue >= Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // T_FT
	// plm.new_transition(events[3], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + zero, Htrue, beta*Htilde, zero}, (counter >= ts) && (Htilde > Hmin) && (Htrue >= Hd) && (Htrue <= Hmax), Ariadne::EventKind::PERMISSIVE); // T_FF

	// If H - deltaH <= Hmin then harvesting is always on, and H will never be > Hmax
	// So every combination starting with F will have two dc at the end F * _ _

	// sign(tau_tl) * abs(newH - Hmin / deltaQ);
	// plm.new_transition(events[4], Ariadne::next({tau_pl, H, H_out, tau_tlc})={abs(mod_tau_tl), Hmin, zero, zero}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (tau_tl >= zero), Ariadne::EventKind::PERMISSIVE); // 
	// plm.new_transition(events[5], Ariadne::next({tau_pl, H, H_out, tau_tlc})={minus_one * abs(mod_tau_tl), Hmin, zero, zero}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (tau_tl < zero), Ariadne::EventKind::PERMISSIVE); // 
	// plm.new_transition(events[5], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FTTF
	// plm.new_transition(events[6], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FTFT
	// plm.new_transition(events[7], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FTFF
	// plm.new_transition(events[8], Ariadne::next({tau_pl, H, H_out, tau_tlc})={zero, Hmin, zero, zero}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ == zero), Ariadne::EventKind::PERMISSIVE); // FFTF
	// plm.new_transition(events[9], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ == zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FFTT
	// plm.new_transition(events[10], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ == zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FFFT
	// plm.new_transition(events[11], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ == zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FFFF

	return plm;
}

Ariadne::HybridAutomaton PLS()
{	
	Ariadne::RealVariable position_slave_d("qs_d");
	Ariadne::RealVariable velocity_slave_d("qs_dot");
	Ariadne::RealVariable position_slave_d_prev("qs_d_prev");
	Ariadne::RealVariable tau_tl("tau_tls");
	Ariadne::RealVariable tau_pl("tau_pls");
	Ariadne::RealVariable deltaH_debug("deltaH_s");

	// Ariadne::RealConstant ts("ts", Ariadne::Decimal(TS));
	// Ariadne::RealVariable counter("cnt");

	Ariadne::RealConstant zero("zero", Ariadne::Decimal(0));
	Ariadne::RealConstant minus_one("minus_one", Ariadne::Decimal(-1.0));
	Ariadne::RealConstant Hmax("Hmax_s", Ariadne::Decimal(HMAX));
	Ariadne::RealConstant Hmin("Hmin_s", Ariadne::Decimal(HMIN));
	Ariadne::RealConstant Hd("Hd_s", Ariadne::Decimal(HD));
	Ariadne::RealConstant beta("beta_s", Ariadne::Decimal(BETA));

	Ariadne::RealVariable H("Hs");
	Ariadne::RealVariable H_in("H+s");
	Ariadne::RealVariable H_out("H-s");

	Ariadne::StringVariable pls_name("pls");
	Ariadne::HybridAutomaton pls(pls_name.name());

	int n_locations = 10;
	int n_events = 9;
	Ariadne::DiscreteLocation loc[n_locations];
	for (size_t i = 0; i < n_locations; i++)
	{
		loc[i] = Ariadne::DiscreteLocation(getPair(pls_name,"s" + std::to_string(i)));
	}
	
	Ariadne::DiscreteEvent clock_event("clock_event");
	Ariadne::DiscreteEvent events[n_events];
	for (size_t i = 0; i < n_events; i++)
	{
		events[i] = Ariadne::DiscreteEvent("pls_t" + std::to_string(i));
	}
	
	Ariadne::RealExpression newH = H + H_in;
	Ariadne::RealExpression deltaQ = (position_slave_d - position_slave_d_prev);
	Ariadne::RealExpression deltaH = deltaQ * tau_tl;
	Ariadne::RealExpression Htilde = newH - deltaH;
	Ariadne::RealExpression mod_tau_tl = newH - Hmin / deltaQ;

	Ariadne::RealExpression Htrue = Htilde - beta * Htilde;
	Ariadne::RealExpression Hfalse = zero;
	 
	for (size_t i = 0; i < n_locations; i++)
	{
		pls.new_mode(loc[i], Ariadne::dot({tau_pl, H, H_out, deltaH_debug}) = {0, 0, 0, 0});
	} 	
	
	pls.new_transition(loc[0], events[0], loc[1], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (Htilde > Hmin), Ariadne::EventKind::PERMISSIVE);
	pls.new_transition(loc[0], events[1], loc[2], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (Htilde <= Hmin), Ariadne::EventKind::PERMISSIVE);

	pls.new_transition(loc[1], events[2], loc[3], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE);
	pls.new_transition(loc[1], events[3], loc[4], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (Htrue <= Hmax), Ariadne::EventKind::PERMISSIVE);

	pls.new_transition(loc[2], events[4], loc[5], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (deltaQ == zero), Ariadne::EventKind::PERMISSIVE);
	pls.new_transition(loc[2], events[5], loc[6], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (deltaQ != zero), Ariadne::EventKind::PERMISSIVE);

	pls.new_transition(loc[6], events[6], loc[7], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (tau_tl > zero), Ariadne::EventKind::PERMISSIVE);
	pls.new_transition(loc[6], events[7], loc[8], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (Ariadne::RealExpression(tau_tl) == zero), Ariadne::EventKind::PERMISSIVE);
	pls.new_transition(loc[6], events[8], loc[9], Ariadne::next({tau_pl, H, H_out, deltaH_debug}) = {tau_pl, H, H_out, deltaH_debug}, (tau_tl < zero), Ariadne::EventKind::PERMISSIVE);

	// Clock transitions
	pls.new_transition(loc[3], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, deltaH_debug})={tau_tl, Hmax, beta*Htilde, deltaH}); // Saturation ON
	pls.new_transition(loc[4], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, deltaH_debug})={tau_tl, Htrue, beta*Htilde, deltaH}); // Saturation OFF
	pls.new_transition(loc[5], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, deltaH_debug})={zero, Hmin, zero, deltaH}); // No energy available and no movement 
	pls.new_transition(loc[7], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, deltaH_debug})={abs(mod_tau_tl), Hmin, zero, deltaH}); // Reducing the torque, sign was positive
	pls.new_transition(loc[8], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, deltaH_debug})={zero, Hmin, zero, deltaH}); // tau_tl was zero, therefore sign is zero
	pls.new_transition(loc[9], clock_event, loc[0], Ariadne::next({tau_pl, H, H_out, deltaH_debug})={minus_one * abs(mod_tau_tl), Hmin, zero, deltaH}); // Reducing the torque, sign was negative

	// pls.new_transition(events[1], Ariadne::next({tau_pl, H, H_out})={tau_tl, Hmax, beta*Htilde}, (counter >= ts) && (Htilde > Hmin) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE);
	// pls.new_transition(events[1], Ariadne::next({tau_pl, H, H_out})={tau_tl, Htrue, beta*Htilde}, (counter >= ts) && (Htilde > Hmin) && (Htrue <= Hmax), Ariadne::EventKind::PERMISSIVE);

	// sign(tau_tl) * abs(newH - Hmin / deltaQ);
	// Ariadne::RealExpression mod_tau_tl = newH - Hmin / deltaQ;
	// pls.new_transition(events[4], Ariadne::next({tau_pl, H, H_out})={abs(mod_tau_tl), Hmin, zero}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (tau_tl >= zero), Ariadne::EventKind::PERMISSIVE); // 
	// pls.new_transition(events[5], Ariadne::next({tau_pl, H, H_out})={minus_one * abs(mod_tau_tl), Hmin, zero}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (tau_tl < zero), Ariadne::EventKind::PERMISSIVE); // 
	// pls.new_transition(events[5], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FTTF
	// pls.new_transition(events[6], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FTFT
	// pls.new_transition(events[7], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ != zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FTFF
	// pls.new_transition(events[8], Ariadne::next({tau_pl, H, H_out})={zero, Hmin, zero}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ == zero), Ariadne::EventKind::PERMISSIVE); // FFTF
	// plm.new_transition(events[9], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ == zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FFTT
	// plm.new_transition(events[10], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ == zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FFFT
	// plm.new_transition(events[11], Ariadne::next({tau_pl, H, H_out, tau_tlc})={tau_tl + tau_tlc_exp, Hmax, beta*Htilde, tau_tlc_exp}, (counter >= ts) && (Htilde <= Hmin) && (deltaQ == zero) && (Htrue < Hd) && (Htrue > Hmax), Ariadne::EventKind::PERMISSIVE); // FFFF

	return pls;
}

Ariadne::HybridAutomaton CommunicationChannel()
{
	// Passing to the other side the discretized values
	Ariadne::RealVariable position_master("qm_d");
	Ariadne::RealVariable velocity_master("qm_dot_d");
	Ariadne::RealVariable position_slave("qs_d");
	Ariadne::RealVariable velocity_slave("qs_dot_d");
	
	Ariadne::RealVariable position_master_m2s("qm_m2s");
	Ariadne::RealVariable velocity_master_m2s("qm_dot_m2s");
	Ariadne::RealVariable position_slave_s2m("qs_s2m");
	Ariadne::RealVariable velocity_slave_s2m("qs_dot_s2m");

	Ariadne::RealVariable H_in_m("H+m");
	Ariadne::RealVariable H_out_m("H-m");
	Ariadne::RealVariable H_in_s("H+s");
	Ariadne::RealVariable H_out_s("H-s");

	Ariadne::HybridAutomaton comm("comm_channel");
	Ariadne::DiscreteLocation loc;

	Ariadne::DiscreteEvent clock_event("clock_event");

	// comm.new_mode(loc, Ariadne::dot({position_master_m2s, velocity_master_m2s, position_slave_s2m, velocity_slave_s2m}) = {0, 0, 0, 0});

	// comm.new_transition(loc, clock_event, loc, Ariadne::next(
	// 	{position_master_m2s, velocity_master_m2s, position_slave_s2m, velocity_slave_s2m}) = 
	// 	{position_master, velocity_master, position_slave, velocity_slave});

	comm.new_mode(loc, Ariadne::let(
		{position_master_m2s, velocity_master_m2s, position_slave_s2m, velocity_slave_s2m, H_in_m, H_in_s}) = 
		{position_master, velocity_master, position_slave, velocity_slave, H_out_s, H_out_m});
	return comm;
}

Ariadne::HybridAutomaton Environment()
{	
	Ariadne::RealConstant position_env("qe", Ariadne::Decimal(QE));
	Ariadne::RealConstant K("K", Ariadne::Decimal(KE));
	Ariadne::RealConstant B("B", Ariadne::Decimal(BE));

	Ariadne::RealVariable position_slave("qs");
	Ariadne::RealVariable velocity_slave("qs_dot");
	Ariadne::RealVariable env_force("fe");

	Ariadne::StringVariable env_name("env");
	Ariadne::HybridAutomaton env(env_name.name());
	Ariadne::DiscreteLocation free_motion(getPair(env_name,"free_motion"));
	Ariadne::DiscreteLocation contact(getPair(env_name,"contact"));

	Ariadne::DiscreteEvent no_force("no_force");
	Ariadne::DiscreteEvent force("force");
	Ariadne::DiscreteEvent no_force_inv("no_force_inv");
	Ariadne::DiscreteEvent force_inv("force_inv");

	env.new_mode(contact, Ariadne::let({env_force}) = {-K * (position_slave - position_env) + B * velocity_slave});
	env.new_mode(free_motion, Ariadne::let({env_force}) = {0});

	env.new_transition(free_motion, force, contact, position_env <= position_slave, Ariadne::EventKind::URGENT);
	env.new_transition(contact, no_force, free_motion, position_env >= position_slave, Ariadne::EventKind::URGENT);
	
	return env;
}

Ariadne::HybridAutomaton TLM()
{
	Ariadne::RealConstant P("P", Ariadne::Decimal(PM));
	Ariadne::RealConstant D("D", Ariadne::Decimal(DM));
	
	Ariadne::RealVariable position_master("qm_d");
	Ariadne::RealVariable velocity_master("qm_dot_d");
	Ariadne::RealVariable position_slave_s2m("qs_s2m");
	Ariadne::RealVariable velocity_slave_s2m("qs_dot_s2m");
	Ariadne::RealVariable torque("tau_tlm");

	Ariadne::HybridAutomaton tlm("tlm");
	Ariadne::DiscreteLocation loc;

	Ariadne::DiscreteEvent clock_event("clock_event");
	
	tlm.new_mode(loc, Ariadne::let({torque}) = {(position_master - position_slave_s2m) * P + (velocity_master - velocity_slave_s2m) * D});

	// tlm.new_mode(loc, Ariadne::dot({torque}) = {0});

	// tlm.new_transition(loc, clock_event, loc, Ariadne::next({torque}) = {(position_master - position_slave_s2m) * P + (velocity_master - velocity_slave_s2m) * D});

	return tlm;
}

Ariadne::HybridAutomaton TLS()
{
	Ariadne::RealConstant P("P", Ariadne::Decimal(PS));
	Ariadne::RealConstant D("D", Ariadne::Decimal(DS));
	
	Ariadne::RealVariable position_master_m2s("qm_m2s");
	Ariadne::RealVariable velocity_master_m2s("qm_dot_m2s");
	Ariadne::RealVariable position_slave("qs_d");
	Ariadne::RealVariable velocity_slave("qs_dot_d");
	Ariadne::RealVariable torque("tau_tls");

	Ariadne::HybridAutomaton tls("tls");
	Ariadne::DiscreteLocation loc;

	Ariadne::DiscreteEvent clock_event("clock_event");

	tls.new_mode(loc, Ariadne::let({torque}) = {(position_master_m2s - position_slave) * P + (velocity_master_m2s - velocity_slave) * D});

	// tls.new_mode(loc, Ariadne::dot({torque}) = {0});
	// tls.new_transition(loc, clock_event, loc, Ariadne::next({torque}) = {(position_master_m2s - position_slave) * P + (velocity_master_m2s - velocity_slave) * D});

	return tls;
}

Ariadne::HybridAutomaton Master()
{
	Ariadne::RealConstant J("J", Ariadne::Decimal(JM));
	Ariadne::RealConstant B("B", Ariadne::Decimal(BM));
	
	Ariadne::RealConstant arm("arm_m", Ariadne::Decimal(ARM_M));
	Ariadne::RealConstant t2c("Kt2c_m", Ariadne::Decimal(KT2CM));
	Ariadne::RealConstant c2v("Kc2v_m", Ariadne::Decimal(KC2VM));

	Ariadne::RealVariable torque("tau_plm");
	Ariadne::RealVariable position_master("qm");
	Ariadne::RealVariable velocity_master("qm_dot");
	Ariadne::RealVariable human_force("h_m_star");

	Ariadne::HybridAutomaton master("master");
	Ariadne::DiscreteLocation loc;

	Ariadne::RealExpression acc_master = (((torque + (human_force * arm)) * t2c * c2v) + (B * velocity_master)) / J; 

	master.new_mode(loc, Ariadne::dot({velocity_master, position_master}) = {acc_master, velocity_master});

	return master;
}

Ariadne::HybridAutomaton Slave()
{
	Ariadne::RealConstant J("J", Ariadne::Decimal(JS));
	Ariadne::RealConstant B("B", Ariadne::Decimal(BS));

	Ariadne::RealConstant arm("arm_s", Ariadne::Decimal(ARM_S));
	Ariadne::RealConstant t2c("Kt2c_s", Ariadne::Decimal(KT2CS));
	Ariadne::RealConstant c2v("Kc2v_s", Ariadne::Decimal(KC2VS));

	Ariadne::RealVariable torque("tau_pls");
	Ariadne::RealVariable env_force("fe");
	Ariadne::RealVariable position_slave("qs");
	Ariadne::RealVariable velocity_slave("qs_dot");

	Ariadne::HybridAutomaton slave("slave");
	Ariadne::DiscreteLocation loc;

	Ariadne::RealExpression acc_slave = (((torque - (env_force * arm)) * t2c * c2v) + (B * velocity_slave)) / J; 

	slave.new_mode(loc, Ariadne::dot({velocity_slave, position_slave}) = {acc_slave, velocity_slave});

	return slave;
}

Ariadne::HybridAutomaton Operator()
{
	Ariadne::RealConstant P("P", Ariadne::Decimal(PO));
	Ariadne::RealConstant D("D", Ariadne::Decimal(DO));
	
	Ariadne::RealVariable position_ref("q_ref");
	Ariadne::RealVariable velocity_ref("q_dot_ref");
	Ariadne::RealVariable position_slave_s2m("qs_s2m");
	Ariadne::RealVariable velocity_slave_s2m("qs_dot_s2m");
	Ariadne::RealVariable human_force("h_m_star");

	Ariadne::HybridAutomaton operator_("human_operator");
	Ariadne::DiscreteLocation loc;

	// Discretized with clock event or continuous?
	operator_.new_mode(loc, Ariadne::let({human_force}) = {(position_ref - position_slave_s2m) * P + (velocity_ref - velocity_slave_s2m) * D});

	return operator_;
}

#if defined (SINE)
Ariadne::HybridAutomaton HumanIntention()
{
	Ariadne::RealConstant freq("freq", Ariadne::Decimal(FREQ));
	Ariadne::RealConstant amp("amp", Ariadne::Decimal(AMP));

	Ariadne::RealVariable position_ref("q_ref");
	Ariadne::RealVariable velocity_ref("q_dot_ref");
	Ariadne::RealVariable t("t");

	Ariadne::HybridAutomaton intention("human_intention");
	Ariadne::DiscreteLocation loc;

	// intention.new_mode(loc, Ariadne::dot({velocity_ref, position_ref,t}) = {1, velocity_ref,1});
	intention.new_mode(loc, Ariadne::dot({velocity_ref, position_ref}) = {amp*2*PI*freq*cos(2*PI*t*freq), velocity_ref});
	// intention.new_mode(loc, Ariadne::let({velocity_ref, position_ref,t}) = {-amp * sin(2*PI*t*freq), velocity_ref,1});
	// intention.new_mode(loc, Ariadne::let({velocity_ref, position_ref}) = {amp * cos(2*PI*t*freq), amp * sin(2*PI*t*freq)}, Ariadne::dot(t) = 1);

	return intention;
}
#elif defined (STEP)
Ariadne::HybridAutomaton HumanIntention()
{
	Ariadne::RealConstant start("start", Ariadne::Decimal(START_TIME));
	Ariadne::RealConstant up("up", Ariadne::Decimal(UP_TIME));
	Ariadne::RealConstant amp("amp", Ariadne::Decimal(AMP));

	Ariadne::RealVariable position_ref("q_ref");
	Ariadne::RealVariable velocity_ref("q_dot_ref");
	Ariadne::RealVariable t("t");

	Ariadne::StringVariable name("human_intention");
	Ariadne::HybridAutomaton intention(name.name());
	Ariadne::DiscreteLocation low(getPair(name,"low"));
	Ariadne::DiscreteLocation changing(getPair(name,"changing"));
	Ariadne::DiscreteLocation stable(getPair(name,"stable"));

	Ariadne::DiscreteEvent change("change");
	Ariadne::DiscreteEvent high("high");

	intention.new_mode(low, Ariadne::dot({velocity_ref, position_ref}) = {0, 0});
	intention.new_mode(changing, Ariadne::dot({velocity_ref, position_ref}) = {0, 0});
	intention.new_mode(stable, Ariadne::dot({velocity_ref, position_ref}) = {0, 0});

	intention.new_transition(low, change, changing, Ariadne::next({position_ref, velocity_ref})={amp, amp}, t >= start, Ariadne::EventKind::URGENT);
	intention.new_transition(changing, high, stable, Ariadne::next({position_ref, velocity_ref})={amp, 0}, t >= up, Ariadne::EventKind::URGENT);

	return intention;
}
#endif

Ariadne::HybridAutomaton Time()
{
	Ariadne::RealVariable t("t");

	Ariadne::HybridAutomaton time("time");
	Ariadne::DiscreteLocation loc;

	time.new_mode(loc, Ariadne::dot({t}) = {1});	
	return time;
}

Ariadne::Void simulate_evolution(const Ariadne::CompositeHybridAutomaton& system, const Nat& log_verbosity)
{
	Ariadne::RealVariable position_master("qm");
	Ariadne::RealVariable velocity_master("qm_dot");
	Ariadne::RealVariable torque_plm("tau_plm");
	Ariadne::RealVariable torque_tlm("tau_tlm");
	Ariadne::RealVariable torque_tlc("tau_tlc");

	Ariadne::RealVariable position_slave("qs");
	Ariadne::RealVariable velocity_slave("qs_dot");
	Ariadne::RealVariable torque_pls("tau_pls");
	Ariadne::RealVariable torque_tls("tau_tls");

	Ariadne::RealVariable position_ref("q_ref");
	Ariadne::RealVariable velocity_ref("q_dot_ref");
	Ariadne::RealVariable human_force("h_m_star");
	
	Ariadne::RealVariable env_force("fe");
	Ariadne::RealVariable t("t");

	Ariadne::RealVariable position_master_d("qm_d");
	Ariadne::RealVariable position_slave_d("qs_d");
	Ariadne::RealVariable velocity_master_d("qm_dot_d");
	Ariadne::RealVariable velocity_slave_d("qs_dot_d");
	Ariadne::RealVariable env_force_d("fe_d");
	Ariadne::RealVariable counter("cnt");

	Ariadne::RealVariable H_m("Hm");
	Ariadne::RealVariable H_in_m("H+m");
	Ariadne::RealVariable H_out_m("H-m");
	Ariadne::RealVariable H_s("Hs");
	Ariadne::RealVariable H_in_s("H+s");
	Ariadne::RealVariable H_out_s("H-s");
	Ariadne::RealVariable deltaHm("deltaH_m");
	Ariadne::RealVariable deltaHs("deltaH_s");

	Ariadne::HybridSimulator simulator;
	simulator.set_step_size(STEP_SIZE);
	simulator.verbosity = log_verbosity;

	Ariadne::HybridRealPoint initial_point(
		{
			getPair("env","free_motion"),
			getPair("human_intention", "low"),
			getPair("pls","s0"),
			getPair("plm","s0")
		},
		{
			position_master = 0, 
			velocity_master = 0,
			position_slave = 0,
			velocity_slave = 0,
			position_ref = 0,
			velocity_ref = 0, 
			t = 0,
			position_master_d = 0,
			position_slave_d = 0,
			velocity_master_d = 0,
			velocity_slave_d = 0,
			env_force_d=0,
			counter = 0,
			torque_plm = 0,
			torque_pls = 0,
			torque_tlc = 0,
			H_m = Ariadne::Decimal( 0.7),
			H_out_m = 0,
			H_s = Ariadne::Decimal( 0.7),
			H_out_s = 0,
			deltaHm = 0,
			deltaHs = 0
		}
	);

	Ariadne::Int tf_c = TF_CONTINUOUS;
	Ariadne::Int tf_d = TF_DISCRETE;
	Ariadne::Real ymin = Ariadne::Decimal(Y_MIN);
	Ariadne::Real ymax = Ariadne::Decimal(Y_MAX);

	Ariadne::HybridTerminationCriterion termination(tf_c, tf_d);

	std::cout << "Computing simulation trajectory...\n" << std::flush;
	auto orbit = simulator.orbit(system,initial_point,termination);
	std::cout << "done!\n" << std::endl;
	
	std::cout << "Logging to file: log.txt" << std::endl;
	std::ofstream log_file;
    log_file.open ("log.txt");
    log_file << orbit << std::endl;
    log_file.close();
	
	std::cout << "Plotting simulation trajectory..\n" << std::flush;
	plot("plots/ref_pos",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=position_ref<=ymax),orbit);
	plot("plots/ref_vel",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=velocity_ref<=ymax),orbit);
	plot("plots/master_pos",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=position_master<=ymax),orbit);
	plot("plots/master_vel",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=velocity_master<=ymax),orbit);
	plot("plots/slave_pos",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=position_slave<=ymax),orbit);
	plot("plots/slave_vel",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=velocity_slave<=ymax),orbit);

	plot("plots/env_force",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=env_force<=ymax),orbit);

	plot("plots/cnt",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, Ariadne::Decimal(0) <=counter<=Ariadne::Decimal(0.005)),orbit);
	plot("plots/qm_d",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=position_master_d<=ymax),orbit);
	plot("plots/qs_d",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=position_slave_d<=ymax),orbit);
	plot("plots/qm_dot_d",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=velocity_master_d<=ymax),orbit);
	plot("plots/qs_dot_d",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=velocity_slave_d<=ymax),orbit);

	plot("plots/torque_pls",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=torque_pls<=ymax),orbit);
	plot("plots/torque_plm",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=torque_plm<=ymax),orbit);
	plot("plots/torque_tls",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=torque_tls<=ymax),orbit);
	plot("plots/torque_tlm",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=torque_tlm<=ymax),orbit);
	plot("plots/torque_tlc",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=torque_tlc<=ymax),orbit);
	plot("plots/H_m",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=H_m<=ymax),orbit);
	plot("plots/H_s",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=H_s<=ymax),orbit);

	plot("plots/H_in_m",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=H_in_m<=ymax),orbit);
	plot("plots/H_out_m",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=H_out_m<=ymax),orbit);
	plot("plots/H_in_s",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=H_in_s<=ymax),orbit);
	plot("plots/H_out_s",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=H_out_s<=ymax),orbit);

	std::cout << "done!\n" << std::endl;
}

Ariadne::Int main(Ariadne::Int argc, const char* argv[])
{	
	load_settings();
	
	Ariadne::Nat log_verbosity = Ariadne::get_verbosity(argc, argv);
	
	Ariadne::CompositeHybridAutomaton system("system", 
	{
		HumanIntention(),
		Operator(),
		Slave(),
		Master(),
		TLM(),
		TLS(),
		Environment(),
		CommunicationChannel(),
		Holder(),
		Clock(),
		Time(),
		OneDelay(),
		PLM(),
		PLS()
	});
	
	Ariadne::CompositeHybridAutomaton::set_default_writer(new Ariadne::CompactCompositeHybridAutomatonWriter());
	std::cout << "\nSystem:\n\n" << system << std::endl;
	simulate_evolution(system,log_verbosity);

	return 0;
}