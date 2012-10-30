/*
PRICE was written by Graham Ruby.
This file is part of PRICE.

PRICE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PRICE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PRICE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ASSEMBLYEXCEPTION_CPP
#define ASSEMBLYEXCEPTION_CPP


#include "AssemblyException.h"

AssemblyException::AssemblyError::AssemblyError(){}
AssemblyException::AssemblyError::AssemblyError(const char msg[]){ _msg = msg; }
const char* AssemblyException::AssemblyError::what() const throw(){ return _msg; }

AssemblyException::ArgError::ArgError(){}
AssemblyException::ArgError::ArgError(const char msg[]){ _msg = msg; }
const char* AssemblyException::ArgError::what() const throw(){ return _msg; }

AssemblyException::FileError::FileError(){}
AssemblyException::FileError::FileError(const char msg[]){ _msg = msg; }
const char* AssemblyException::FileError::what() const throw(){ return _msg; }

AssemblyException::LogicError::LogicError(){}
AssemblyException::LogicError::LogicError(const char msg[]){ _msg = msg; }
const char* AssemblyException::LogicError::what() const throw(){ return _msg; }

AssemblyException::LengthError::LengthError(){}
AssemblyException::LengthError::LengthError(const char msg[]){ _msg = msg; }
const char* AssemblyException::LengthError::what() const throw(){ return _msg; }

AssemblyException::ImplementationError::ImplementationError(){}
AssemblyException::ImplementationError::ImplementationError(const char msg[]){ _msg = msg; }
const char* AssemblyException::ImplementationError::what() const throw(){ return _msg; }

AssemblyException::CallingError::CallingError(){}
AssemblyException::CallingError::CallingError(const char msg[]){ _msg = msg; }
const char* AssemblyException::CallingError::what() const throw(){ return _msg; }

AssemblyException::StateError::StateError(){}
AssemblyException::StateError::StateError(const char msg[]){ _msg = msg; }
const char* AssemblyException::StateError::what() const throw(){ return _msg; }

AssemblyException::CopyConstructorError::CopyConstructorError(){}
AssemblyException::CopyConstructorError::CopyConstructorError(const char msg[]){ _msg = msg; }
const char* AssemblyException::CopyConstructorError::what() const throw(){ return _msg; }

AssemblyException::CopyAssignmentOperatorError::CopyAssignmentOperatorError(){}
AssemblyException::CopyAssignmentOperatorError::CopyAssignmentOperatorError(const char msg[]){ _msg = msg; }
const char* AssemblyException::CopyAssignmentOperatorError::what() const throw(){ return _msg; }

AssemblyException::KmerN_Error::KmerN_Error(){}
AssemblyException::KmerN_Error::KmerN_Error(const char msg[]){ _msg = msg; }
const char* AssemblyException::KmerN_Error::what() const throw(){ return _msg; }

#endif
