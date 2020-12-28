#include <iostream>
#include <ariadne.hpp>

#define PROTOTYPE 1

template<typename T>
void print(T in, bool newline = true)
{
	if (newline)
		std::cout << in << std::endl;
	else
		std::cout << in;
}

#if PROTOTYPE
	#include "yaml-cpp/yaml.h"

	std::string config_path = "../config.yml";

	#define PI Ariadne::pi
	// ENVIRONMENT
	double QE, KE, BE;
	// OPERATOR
	double PO, DO, FREQ, AMP;
	// SLAVE
	double ARM_S, JS, BS, BhS, JhS, KT2CS, KC2VS, PS, DS;
	// MASTER
	double ARM_M, JM, BM, BhM, JhM, KT2CM, KC2VM, PM, DM;
	// SIMULATION
	Ariadne::Int TF_CONTINUOUS, TF_DISCRETE;
	// PLOTS
	double Y_MIN, Y_MAX;

	void load_settings()
	{
		YAML::Node config = YAML::LoadFile(config_path);
		QE = config["QE"].as<double>(); 
		KE = config["KE"].as<double>(); 
		BE = config["BE"].as<double>(); 
		PO = config["PO"].as<double>();
		DO = config["DO"].as<double>(); 
		FREQ = config["FREQ"].as<double>(); 
		AMP = config["AMP"].as<double>(); 	
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

		if (config["DEBUG"].as<bool>()) {
			print("ENVIRONMENT");
			print("QE: " + std::to_string(QE));
			print("KE: " + std::to_string(KE));
			print("BE: " + std::to_string(BE));
			print("\nOPERATOR");
			print("PO: " + std::to_string(PO));
			print("DO: " + std::to_string(DO));
			print("FREQ: " + std::to_string(FREQ));
			print("AMP: " + std::to_string(AMP));
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
			print("TF_CONTINUOUS: " + std::to_string(TF_CONTINUOUS));
			print("TF_DISCRETE: " + std::to_string(TF_DISCRETE));
			print("\nPLOTS");
			print("Y_MIN: " + std::to_string(Y_MIN));
			print("Y_MAX: " + std::to_string(Y_MAX));
			print("");
		}
	}
#else
	#include "defines.hpp"
#endif

// Mimic the overriden operator outside the namespace 
Ariadne::Pair<Ariadne::StringVariable, Ariadne::String> getPair(Ariadne::StringVariable base, std::string location) {{ using namespace Ariadne; return base|location;};}
Ariadne::Pair<Ariadne::StringVariable, Ariadne::String> getPair(std::string base, std::string location) {{ using namespace Ariadne; return Ariadne::StringVariable("base")|Ariadne::StringConstant(location);};}

Ariadne::HybridAutomaton CommunicationChannel()
{
	Ariadne::RealVariable position_master("qm");
	Ariadne::RealVariable velocity_master("qm_dot");
	Ariadne::RealVariable position_slave("qs");
	Ariadne::RealVariable velocity_slave("qs_dot");
	
	Ariadne::RealVariable position_master_m2s("qm_m2s");
	Ariadne::RealVariable velocity_master_m2s("qm_dot_m2s");
	Ariadne::RealVariable position_slave_s2m("qs_s2m");
	Ariadne::RealVariable velocity_slave_s2m("qs_dot_s2m");

	Ariadne::HybridAutomaton comm("comm_channel");
	Ariadne::DiscreteLocation loc;

	comm.new_mode(loc, Ariadne::let(
		{position_master_m2s, velocity_master_m2s, position_slave_s2m, velocity_slave_s2m}) = 
		{position_master, velocity_master, position_slave, velocity_slave}
		);

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

	env.new_mode(free_motion, Ariadne::let({env_force}) = {K * (position_slave - position_env) + B * velocity_slave});
	env.new_mode(contact, Ariadne::let({env_force}) = {0});

	Ariadne::RealConstant delta("delta", Ariadne::Decimal(0.01));

	// TODO ?
	// env.new_invariant(free_motion, position_slave >= position_env, force_inv);
	// env.new_invariant(contact, position_env < position_slave, no_force_inv);

	env.new_transition(free_motion, force, contact, position_env <= position_slave + delta, Ariadne::EventKind::URGENT);
	env.new_transition(contact, no_force, free_motion, position_env >= position_slave - delta, Ariadne::EventKind::URGENT);
	
	return env;
}

Ariadne::HybridAutomaton TLM()
{
	Ariadne::RealConstant P("P", Ariadne::Decimal(PM));
	Ariadne::RealConstant D("D", Ariadne::Decimal(DM));
	
	Ariadne::RealVariable position_master("qm");
	Ariadne::RealVariable velocity_master("qm_dot");
	Ariadne::RealVariable position_slave_s2m("qs_s2m");
	Ariadne::RealVariable velocity_slave_s2m("qs_dot_s2m");
	Ariadne::RealVariable torque("tau_m");

	Ariadne::HybridAutomaton tlm("tlm");
	Ariadne::DiscreteLocation loc;

	tlm.new_mode(loc, Ariadne::let({torque}) = {(position_master - position_slave_s2m) * P + (velocity_master - velocity_slave_s2m) * D});

	return tlm;
}

Ariadne::HybridAutomaton TLS()
{
	Ariadne::RealConstant P("P", Ariadne::Decimal(PS));
	Ariadne::RealConstant D("D", Ariadne::Decimal(DS));
	
	Ariadne::RealVariable position_master_m2s("qm_m2s");
	Ariadne::RealVariable velocity_master_m2s("qm_dot_m2s");
	Ariadne::RealVariable position_slave("qs");
	Ariadne::RealVariable velocity_slave("qs_dot");
	Ariadne::RealVariable torque("tau_s");

	Ariadne::HybridAutomaton tls("tls");

	Ariadne::DiscreteLocation loc;

	tls.new_mode(loc, Ariadne::let({torque}) = {(position_master_m2s - position_slave) * P + (velocity_master_m2s - velocity_slave) * D});

	return tls;
}

Ariadne::HybridAutomaton Master()
{
	Ariadne::RealConstant J("J", Ariadne::Decimal(JM));
	Ariadne::RealConstant B("B", Ariadne::Decimal(BM));
	
	Ariadne::RealConstant arm("arm_m", Ariadne::Decimal(ARM_M));
	Ariadne::RealConstant t2c("Kt2c_m", Ariadne::Decimal(KT2CM));
	Ariadne::RealConstant c2v("Kc2v_m", Ariadne::Decimal(KC2VM));

	Ariadne::RealVariable torque("tau_m");
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

	Ariadne::RealVariable torque("tau_s");
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
	
	Ariadne::RealVariable position_ref("qm_ref");
	Ariadne::RealVariable velocity_ref("qm_dot_ref");
	Ariadne::RealVariable position_slave_s2m("qs_s2m");
	Ariadne::RealVariable velocity_slave_s2m("qs_dot_s2m");
	Ariadne::RealVariable human_force("h_m_star");

	Ariadne::HybridAutomaton operator_("human_operator");
	Ariadne::DiscreteLocation loc;

	operator_.new_mode(loc, Ariadne::let({human_force}) = {(position_ref - position_slave_s2m) * P + (velocity_ref - velocity_slave_s2m) * D});

	return operator_;
}

Ariadne::HybridAutomaton HumanIntention()
{
	Ariadne::RealConstant freq("freq", Ariadne::Decimal(FREQ));
	Ariadne::RealConstant amp("amp", Ariadne::Decimal(AMP));
	
	Ariadne::RealVariable position_ref("qm_ref");
	Ariadne::RealVariable velocity_ref("qm_dot_ref");
	Ariadne::RealVariable t("t");

	Ariadne::HybridAutomaton intention("human_intention");
	Ariadne::DiscreteLocation loc;
	
	// intention.new_mode(loc, Ariadne::dot({velocity_ref, position_ref,t}) = {1, velocity_ref,1});
	intention.new_mode(loc, Ariadne::dot({velocity_ref, position_ref, t}) = {amp*2*PI*freq*cos(2*PI*t*freq), velocity_ref, 1});
	// intention.new_mode(loc, Ariadne::let({velocity_ref, position_ref,t}) = {-amp * sin(2*PI*t*freq), velocity_ref,1});
	// intention.new_mode(loc, Ariadne::let({velocity_ref, position_ref}) = {amp * cos(2*PI*t*freq), amp * sin(2*PI*t*freq)}, Ariadne::dot(t) = 1);
	
	return intention;
}

Ariadne::Void simulate_evolution(const Ariadne::CompositeHybridAutomaton& system, const Nat& log_verbosity)
{
	Ariadne::RealVariable position_master("qm");
	Ariadne::RealVariable velocity_master("qm_dot");
	Ariadne::RealVariable torque_m("tau_m");

	Ariadne::RealVariable position_slave("qs");
	Ariadne::RealVariable velocity_slave("qs_dot");
	Ariadne::RealVariable torque_s("tau_s");
	
	Ariadne::RealVariable position_master_m2s("qm_m2s");
	Ariadne::RealVariable velocity_master_m2s("qm_dot_m2s");
	Ariadne::RealVariable position_slave_s2m("qs_s2m");
	Ariadne::RealVariable velocity_slave_s2m("qs_dot_s2m");

	Ariadne::RealVariable position_ref("qm_ref");
	Ariadne::RealVariable velocity_ref("qm_dot_ref");
	Ariadne::RealVariable human_force("h_m_star");
	
	Ariadne::RealVariable env_force("fe");
	Ariadne::RealVariable t("t");

	Ariadne::HybridSimulator simulator;
	simulator.set_step_size(0.01);
	simulator.verbosity = log_verbosity;

	Ariadne::HybridRealPoint initial_point(
		{
			getPair("env","free_motion")
		},
		{
			position_master=0, 
			velocity_master=0,
			position_slave=0,
			velocity_slave=0,
			position_ref=0,
			velocity_ref=0, 
			t = 0
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

	std::cout << "Plotting simulation trajectory..\n" << std::flush;
	plot("ref_pos",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=position_ref<=ymax),orbit);
	plot("ref_vel",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=velocity_ref<=ymax),orbit);
	plot("master_pos",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=position_master<=ymax),orbit);
	plot("master_vel",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=velocity_master<=ymax),orbit);
	plot("slave_pos",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=position_slave<=ymax),orbit);
	plot("slave_vel",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=velocity_slave<=ymax),orbit);

	plot("torque_m",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=torque_m<=ymax),orbit);
	plot("torque_s",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=torque_s<=ymax),orbit);
	plot("env_force",Ariadne::Axes2d(0<=Ariadne::TimeVariable()<= tf_c, ymin <=env_force<=ymax),orbit);

	std::cout << "done!\n" << std::endl;
}

Ariadne::Int main(Ariadne::Int argc, const char* argv[])
{	
	#if PROTOTYPE
		print("Loading parameters from yaml");
		load_settings();
	#else
		print("Loading parameters from header");
	#endif

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
		CommunicationChannel()
	});
	
	Ariadne::CompositeHybridAutomaton::set_default_writer(new Ariadne::CompactCompositeHybridAutomatonWriter());
	std::cout << "\nSystem:\n\n" << system << std::endl;
	simulate_evolution(system,log_verbosity);

	return 0;
}