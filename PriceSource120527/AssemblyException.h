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

/* Exception classes that are specific to the Assembler.  There will be
one parent exception class (AssemblyError) that will encompass all 
Assembler-derived exceptions (allowing us to distinguish between those
exceptions that we chose to throw and specified and therefore might want 
to catch versus those thrown by C++ STL components that reflect an
unintentional programming error on our part). */

#ifndef ASSEMBLYEXCEPTION_H
#define ASSEMBLYEXCEPTION_H

#include <exception>
using namespace std;


namespace AssemblyException {

  /* The system-specific parent exception. */
  class AssemblyError : public exception{
  public:
    AssemblyError();
    AssemblyError(const char msg[]);
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };


  /* Bad arg into a function */
  class ArgError : public AssemblyError{
  public:
    ArgError();
    ArgError(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };

  /* File trouble */
  class FileError : public AssemblyError{
  public:
    FileError();
    FileError(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };

  /* Some sort of logical inconsistancy */
  class LogicError : public AssemblyError{
  public:
    LogicError();
    LogicError(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };

  /* Something is the wrong length */
  class LengthError : public AssemblyError{
  public:
    LengthError();
    LengthError(const char msg[]);
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };

  /* Something has not been implemented */
  class ImplementationError : public AssemblyError{
  public:
    ImplementationError();
    ImplementationError(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };

  /* Something has been called that shouldn't be (a pre-function-call check was skipped) */
  class CallingError : public AssemblyError{
  public:
    CallingError();
    CallingError(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };

  /* invalid internal state of some object */
  class StateError : public AssemblyError{
  public:
    StateError();
    StateError(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };

  /* The copy constructor has been called on an object that should only
   * be passed by reference.
   */
  class CopyConstructorError : public AssemblyError{
  public:
    CopyConstructorError();
    CopyConstructorError(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };

  /* The copy assignment operator has been called on an object that should only
   * be passed by reference.
   */
  class CopyAssignmentOperatorError : public AssemblyError{
  public:
    CopyAssignmentOperatorError();
    CopyAssignmentOperatorError(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };


  /* a kmer string with an N in it */
  class KmerN_Error : public ArgError{
  public:
    KmerN_Error();
    KmerN_Error(const char msg[]);
    char* msg();
    virtual const char* what() const throw();
  private:
    const char* _msg;
  };



}


#endif
