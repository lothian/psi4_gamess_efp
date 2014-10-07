#ifndef GAMESSPARSER_H_
#define GAMESSPARSER_H_

#include <iosfwd>
#include "libmints/typedefs.h"

namespace psi{

class Options;

class GamessOutputParser{
  /**
   * Parses a GAMESS output for core Hamiltonian and MO coefficients, ordering
   * them and renormalizing them for a Psi4 calculation.
   */
  protected:
    /// The options object
    Options &options_;
    /// The C matrix as parsed from GAMESS
    SharedMatrix CGamess_;
    /// The C matrix, ordered and normalized for Psi4
    SharedMatrix C_;
    /// The GAMESS->Psi4 AO transformation matrix for the H integrals
    SharedMatrix UH_;
    /// The GAMESS->Psi4 AO transformation matrix for the MO coefficients
    SharedMatrix UC_;
    /// The H matrix, as parsed from GAMESS
    SharedMatrix HGamess_;
    /// The H matrix, ordered and normalized for Psi4
    SharedMatrix H_;
    /// The number of MOs found in the MO coefficients
    int nmo_;
    /// The number of SOs found in the H and C matrices
    int nso_;

    /// H matrix included in the output
    bool H_found_; 

    /// Reads a stream, to parse for MO coefficients
    void parse_mos(std::ifstream &gamessout);
    /// Reads a stream, to parse for the core Hamiltonian
    void parse_H(std::ifstream &gamessout);
    /// Forms the GAMESS->Psi4 reordering / renormalization matrix
    void build_U_and_rotate();
  public:
    GamessOutputParser(Options &options);
    /// Get the C matrix from GAMESS
    SharedMatrix C() const { return C_; }
    /// Get the H matrix from GAMESS
    SharedMatrix H() const { return H_; }
};

} // End namespace

#endif
