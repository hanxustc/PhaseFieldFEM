#include "nonlocalApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
nonlocalApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

nonlocalApp::nonlocalApp(const InputParameters & parameters) : MooseApp(parameters)
{
  nonlocalApp::registerAll(_factory, _action_factory, _syntax);
}

nonlocalApp::~nonlocalApp() {}

void
nonlocalApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<nonlocalApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"nonlocalApp"});
  Registry::registerActionsTo(af, {"nonlocalApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
nonlocalApp::registerApps()
{
  registerApp(nonlocalApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
nonlocalApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  nonlocalApp::registerAll(f, af, s);
}
extern "C" void
nonlocalApp__registerApps()
{
  nonlocalApp::registerApps();
}
