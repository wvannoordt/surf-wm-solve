#include "HyWall.h"
#include "geolytical.h"
#include "PTL.h"
#include "print.h"

#include <sstream>
struct input_data_t
{
	double dist;
};
bool streq(const std::string& s1, const std::string& s2) {return s1.compare(s2)==0;}
void ReadInput(const std::string& filename, decltype(HyWall::settings)& settings, input_data_t& in_data, bool& ran_init)
{
	PTL::PropertyTree input;
	ran_init = false;
	bool failed = false;
	bool is_init = false;
	std::string failed_message;
	try
	{
		input.Read(filename);
	}
	catch (PTL::PTLException e)
	{
		failed = true;
		is_init = filename == "--init";
		if (!is_init)
		{
			failed_message = e.what();
			std::cout << "FAILURE MESSAGE:\n" << failed_message << std::endl;
		}
	}
	
	settings.enableWallModel = true;
	settings.readRestart = false;
	input["WallModel"]["rayDim"].MapTo(&settings.rayDim)                   = new PTL::PTLInteger(30, "number of ray points");
	settings.asyncSolve = false;
	input["WallModel"]["verboseLevel"].MapTo(&settings.verboseLevel)          = new PTL::PTLInteger(1, "debug output level");
	input["WallModel"]["maxIterations"].MapTo(&settings.maxIterations)        = new PTL::PTLInteger(100, "Max. iterations");
	input["WallModel"]["wallSpacing"].MapTo(&settings.wallSpacing)            = new PTL::PTLDouble(1e-6, "Max. iterations");
	input["WallModel"]["wallTemperature"].MapTo(&settings.wallTemperature)    = new PTL::PTLDouble(100, "Wall Temperature");
	input["WallModel"]["adiabaticWall"].MapTo(&settings.adiabaticWall)        = new PTL::PTLBoolean(true, "Adiabatic wall");
	input["WallModel"]["fluidCp"].MapTo(&settings.fluidCp)                    = new PTL::PTLDouble(1005.0, "Specific heat");
	input["WallModel"]["turbPradntl"].MapTo(&settings.turbPradntl)            = new PTL::PTLDouble(0.72, "Turbulent Prandtl");
	input["WallModel"]["fluidPrandtl"].MapTo(&settings.fluidPrandtl)          = new PTL::PTLDouble(0.9, "Laminar Prandtl");
	input["WallModel"]["vanDriestAPlus"].MapTo(&settings.vanDriestAPlus)      = new PTL::PTLDouble(17.0, "van Driest Constant");
	input["WallModel"]["gasConstant"].MapTo(&settings.gasConstant)            = new PTL::PTLDouble(287.0, "Gas constant");
	input["WallModel"]["enableTransitionSensor"].MapTo(&settings.enableTransitionSensor) = new PTL::PTLBoolean(false, "Enable Transition Sensor");
	
	std::string mom_eq_str, trb_eq_str, eng_eq_str;
	input["WallModel"]["momentumEquationType"].MapTo(&mom_eq_str)          = new PTL::PTLString("ODE", "Momentum equation type");
	input["WallModel"]["turbulenceEquationType"].MapTo(&trb_eq_str)        = new PTL::PTLString("vanDriest", "Turbulence equation type");
	input["WallModel"]["energyEquationType"].MapTo(&eng_eq_str)            = new PTL::PTLString("ODE", "Energy equation type");
	
	input["WallModel"]["momentumUnderRelaxationODE"].MapTo(&settings.momentumUnderRelaxationODE)     = new PTL::PTLDouble(0.8, "Relaxation factor for momentum ODE");
	input["WallModel"]["turbulenceUnderRelaxationODE"].MapTo(&settings.turbulenceUnderRelaxationODE) = new PTL::PTLDouble(0.6, "Relaxation factor for turbulence ODE");
	input["WallModel"]["energyUnderRelaxationODE"].MapTo(&settings.energyUnderRelaxationODE)         = new PTL::PTLDouble(0.7, "Relaxation factor for energy ODE");
	input["WallModel"]["includeMomentumRhs"].MapTo(&settings.includeMomentumRhs)                     = new PTL::PTLBoolean(false, "Include the parameterized convection term");
	input["WallModel"]["isCompressible"].MapTo(&settings.isCompressible)                             = new PTL::PTLBoolean(false, "Use variable density");
	input["WallModel"]["suthViscRef"].MapTo(&settings.suthViscRef)                                   = new PTL::PTLDouble(1.45151376745308e-06, "Reference viscosity for viscosity power law");
	input["WallModel"]["suthTRef"].MapTo(&settings.suthTRef)                                         = new PTL::PTLDouble(110.4, "Reference temperature for viscosity power law");
	
	std::string visc_law_str;
	input["WallModel"]["viscousLaw"].MapTo(&visc_law_str)                                            = new PTL::PTLString("constant", "Viscous law");
	
	if (!failed) input.StrictParse();
	     if (streq(mom_eq_str, "allmaras")) {settings.momentumEquationType = HyCore::momentum::allmaras;}
	else if (streq(mom_eq_str, "ODE"))      {settings.momentumEquationType = HyCore::momentum::ODE;}
	else if (!failed) {std::cout << "Invalid value for momentum equation" << std::endl; abort();}

	     if (streq(trb_eq_str, "linear"))    {settings.turbulenceEquationType = HyCore::turbulence::linear;}
	else if (streq(trb_eq_str, "ODE"))       {settings.turbulenceEquationType = HyCore::turbulence::ODE;}
	else if (streq(trb_eq_str, "vanDriest")) {settings.turbulenceEquationType = HyCore::turbulence::vanDriest;}
	else if (!failed) {std::cout << "Invalid value for turbulence equation" << std::endl; abort();}

	     if (streq(eng_eq_str, "croccoBusemann")) {settings.energyEquationType = HyCore::energy::croccoBusemann;}
	else if (streq(eng_eq_str, "ODE"))            {settings.energyEquationType = HyCore::energy::ODE; }
	else if (streq(eng_eq_str, "linear"))         {settings.energyEquationType = HyCore::energy::linear;}
	else if (!failed) {std::cout << "Invalid value for energy equation" << std::endl; abort();}

	     if (streq(visc_law_str, "constant"))   {settings.viscousLaw = HyCore::visclaw::constant;}
	else if (streq(visc_law_str, "sutherland")) {settings.viscousLaw = HyCore::visclaw::sutherland;}
	else if (streq(visc_law_str, "PowerLaw"))   {settings.viscousLaw = HyCore::visclaw::PowerLaw;}
	else if (!failed) {std::cout << "Invalid value for energy equation" << std::endl; abort();}
	else if (is_init)
	{
		ran_init = true;
		std::cout << "Initializing input file..." << std::endl;
		input.CreateDefaultValuesFile("input.ptl");
		std::cout << "See \"input.ptl\"" << std::endl;
		exit(0);
	}
	else
	{
		std::cout << "Cannot find input file \"" << filename << "\". Try running with first argument \"--init\" instead." << std::endl;
		abort();
	}
}

struct arg_t
{
    arg_t(const std::string& raw_in) { raw = raw_in; }
    std::string raw;
    template <typename out_t> operator out_t() const
    {
        out_t output;
        std::stringstream ss;
        ss << raw;
        ss >> output;
        return output;
    }
};

std::string get_string(const char* cc)
{
    return std::string(cc);
}

template <typename data_t> std::string get_string(const data_t& t)
{
    return std::to_string(t);
}

template <typename def_t> arg_t get_input_arg(int argc, char** argv, const int& position, const def_t& default_initializer)
{
    std::string raw_default = get_string(default_initializer);
    if (argc+1 > position)
    {
        return arg_t(argv[position+1]);
    }
    else
    {
        return arg_t(raw_default);
    }
}

int main(int argc, char** argv)
{
    std::string ptl_file = get_input_arg(argc, argv, 0, "input.ptl");
	MPI_Init(&argc, &argv);
	bool was_init = false;
	HyWall::Initialize(MPI_COMM_WORLD, 4);
    input_data_t input_data;
	ReadInput(ptl_file, HyWall::settings, input_data, was_init);
	HyWall::SetDomainSize(1);
	HyWall::DefineVariables();
	std::vector<double> ff_in(6);
	ff_in[0] = input_data.p;
	ff_in[1] = input_data.u;
	ff_in[2] = input_data.v;
	ff_in[3] = input_data.w;
	ff_in[4] = input_data.T;
	ff_in[5] = input_data.mu_t;
	HyWall::PassFlowfieldVariables(&ff_in[0], 1);
	HyWall::PassVariable("in:distance",    &input_data.dist);
	HyWall::PassVariable("in:x",           &input_data.x);
	HyWall::PassVariable("in:y",           &input_data.y);
	HyWall::PassVariable("in:z",           &input_data.z);
	HyWall::PassVariable("in:rho",         &input_data.rho);
	HyWall::PassVariable("in:mu_lam",      &input_data.mu);
	HyWall::PassVariable("in:momRHS",      &input_data.mom_rhs);
	HyWall::PassVariable("in:dpdx",        &input_data.dp_dx);
	double dummy = 1.0;
	if (HyWall::settings.enableTransitionSensor)
	{
		HyWall::SetTimeStep(1.0);
		HyWall::PassVariable("aux:strain_rate",    &dummy);
		HyWall::PassVariable("aux:sensor_preMult", &dummy);
	}
    
	HyWall::PassVariable("out:vorticity",    &output_data.vort);
	HyWall::PassVariable("out:tau",          &output_data.tau);
	HyWall::PassVariable("out:heatflux",     &output_data.qw);
	HyWall::PassVariable("out:failurelevel", &output_data.fail);
	HyWall::Allocate();
	HyWall::Solve();
	HyWall::Finalize();
	MPI_Finalize();
	return 0;
}