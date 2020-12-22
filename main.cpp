#include <iostream>

#include <ariadne.hpp>

#include "defines.hpp"

// file:///home/kriato/ariadne/build/docs/html/index.html

template<typename T>
void print(T in, bool newline = true)
{
	if (newline)
		std::cout << in << std::endl;
	else
		std::cout << in;
}



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

	Ariadne::HybridAutomaton env("env");
	Ariadne::DiscreteLocation loc;

	env.new_mode(loc, Ariadne::let({env_force}) = {K * (position_slave - position_env) + B * velocity_slave});

	return env;
}

Ariadne::HybridAutomaton TLM()
{
	Ariadne::RealConstant P("P", PM);
	Ariadne::RealConstant D("D", DM);
	
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
	Ariadne::RealConstant P("P", PS);
	Ariadne::RealConstant D("D", DS);
	
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
	Ariadne::RealConstant Bh("Bh", Ariadne::Decimal(BhM));
	Ariadne::RealConstant Jh("Jh", Ariadne::Decimal(JhM));
	// is this the same as _decimal?
	Ariadne::RealConstant arm("arm_m", Ariadne::Decimal(ARM_M));
	Ariadne::RealConstant t2c("Kt2c_m", Ariadne::Decimal(KT2CM));
	Ariadne::RealConstant c2v("Kc2v_m", Ariadne::Decimal(KC2VM));

	Ariadne::RealVariable torque("tau_m");
	Ariadne::RealVariable position_master("qm");
	Ariadne::RealVariable velocity_master("qm_dot");
	Ariadne::RealVariable human_force("h_m_star");

	Ariadne::HybridAutomaton master("master");
	Ariadne::DiscreteLocation loc;
																												  // accelerazione
	master.new_mode(loc, Ariadne::dot({velocity_master, position_master}) = {((torque - ((velocity_master * Bh + velocity_master * Jh) - human_force) * arm) * t2c * c2v) / (J), velocity_master});

	return master;
}

Ariadne::HybridAutomaton Slave()
{
	Ariadne::RealConstant J("J", Ariadne::Decimal(JS));
	Ariadne::RealConstant B("B", Ariadne::Decimal(BS));
	// is this the same as _decimal?
	Ariadne::RealConstant arm("arm_s", Ariadne::Decimal(ARM_S));
	Ariadne::RealConstant t2c("Kt2c_s", Ariadne::Decimal(KT2CS));
	Ariadne::RealConstant c2v("Kc2v_s", Ariadne::Decimal(KC2VS));

	Ariadne::RealVariable torque("tau_s");
	Ariadne::RealVariable env_force("fe");
	Ariadne::RealVariable position_slave("qs");
	Ariadne::RealVariable velocity_slave("qs_dot");

	Ariadne::HybridAutomaton slave("slave");
	Ariadne::DiscreteLocation loc;

	slave.new_mode(loc, Ariadne::dot({velocity_slave, position_slave}) = {((torque - (env_force * arm)) * t2c * c2v) + B * velocity_slave / (J), velocity_slave});

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
	Ariadne::RealConstant freq("freq", FREQ);
	Ariadne::RealConstant amp("amp", AMP);
	
	Ariadne::RealVariable position_ref("qm_ref");
	Ariadne::RealVariable velocity_ref("qm_dot_ref");
	Ariadne::RealVariable t("t");

	Ariadne::HybridAutomaton intention("human_intention");
	Ariadne::DiscreteLocation loc;

	//intention.new_mode(loc, Ariadne::dot({velocity_ref, position_ref,t}) = {1, velocity_ref,1});
	// intention.new_mode(loc, Ariadne::dot({velocity_ref, position_ref,t}) = {-amp * sin(2*PI*t*freq), amp * cos(2*PI*t*freq),1});
	// intention.new_mode(loc, Ariadne::let({velocity_ref, position_ref,t}) = {-amp * sin(2*PI*t*freq), velocity_ref,1});
	intention.new_mode(loc, Ariadne::let({velocity_ref, position_ref}) = {amp * cos(2*PI*t*freq), amp * sin(2*PI*t*freq)}, Ariadne::dot(t) = 1);

	// intention.new_mode(loc, Ariadne::dot({t}) = {1});
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

	Ariadne::HybridSimulator simulator;
	simulator.set_step_size(0.01);
	simulator.verbosity = log_verbosity;

	Ariadne::HybridRealPoint initial_point({},
	{
		position_master=0, 
		velocity_master=0,
		torque_m=0, 
		position_slave=0,
		velocity_slave=0, 
		torque_s=0,
		position_ref=0, 
		velocity_ref=0,
		human_force=0,
		env_force=0
	});

	Ariadne::Int tf_c = TF_CONTINUOUS;
	Ariadne::Int tf_d = TF_DISCRETE;
	Ariadne::Real ymin = Y_MIN;
	Ariadne::Real ymax = Y_MAX;

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
	std::cout << "done!\n" << std::endl;
}

Ariadne::Int main(Ariadne::Int argc, const char* argv[])
{	
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
	std::cout << "System:\n" << system << std::endl;
	simulate_evolution(system,log_verbosity);

	return 0;
}