/** Exceptions thrown from TrajectoryParametrization and dependent subsystems.
 */

#ifndef Propagation_Exceptions_H
#define Propagation_Exceptions_H

//#include <exception>
#include "FWCore/Utilities/interface/Exception.h"

#include <string>

/// Common base class

class PropagationException : public cms::Exception {
public:
  PropagationException( const std::string& message) throw() :  cms::Exception(message)  {}
  virtual ~PropagationException() throw() {}
private:
};

#endif
