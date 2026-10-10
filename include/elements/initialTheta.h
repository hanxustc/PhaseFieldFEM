#pragma once
#include "NodalUserObject.h"
#include <vector>

class initialTheta : public NodalUserObject
{
public:
  static InputParameters validParams();
  initialTheta(const InputParameters &);

  const std::vector<Real> & getThetaInit() const;
  std::vector<Real> _theta_init;

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & other) override;
};
