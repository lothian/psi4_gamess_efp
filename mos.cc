#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

INIT_PLUGIN

using namespace boost;
using namespace std;

namespace psi{ namespace psi4_gamess_mos {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "PSI4_GAMESS_MOS"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" 
PsiReturnType psi4_gamess_mos(Options& options)
{
    int print = options.get_int("PRINT");

    boost::shared_ptr<PSIO> psio(_default_psio_lib_);
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");
    if(wfn->nirrep() != 1) throw PSIEXCEPTION("Must be run w/o symmetry for now.");
    boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));
    fprintf(outfile, "Nuclear repulsion energy = %20.12f\n", chkpt->rd_enuc());
    int nao, nso, nmo;
    fprintf(outfile, "Number of AOs = %d\n", nao = chkpt->rd_nao());
    fprintf(outfile, "Number of SOs = %d\n", nso = chkpt->rd_nso());
    fprintf(outfile, "Number of MOs = %d\n", nmo = chkpt->rd_nmo());

    SharedMatrix scf_SO = wfn->Ca_subset("SO", "ALL"); // AO vs. SO definition in PSI4 is screwed up
    // scf_SO->print(outfile);

    boost::shared_ptr<BasisSet> basisset = wfn->basisset();
    boost::shared_ptr<IntegralFactory> integral = wfn->integral();
    PetiteList petite(basisset, integral, true);
    SharedMatrix aotoso = petite.aotoso();
    // aotoso->print(outfile);

    SharedMatrix scf_AO(new Matrix("SCF AO x MO", nao, nmo));
    scf_AO->gemm(false, false, 1.0, aotoso, scf_SO, 0.0);
    // scf_AO->print(outfile);   

    // Read MO data from GAMESS output file
    ifstream input;
    input.open("mos.txt");
    if(!input.good()) throw PSIEXCEPTION("Error opening mos.txt file.");
    char buf[512];
    SharedVector energies = SharedVector(new Vector("MO Energies", nmo));
    while(!input.eof()) {
      input.getline(buf, 512);  // grab a line from input
      stringstream cppbuf(buf);
      string token;
      vector<string> tokens;
      while(cppbuf >> token) tokens.push_back(token);
      if(tokens.size()) { // skip blank lines
        
        fprintf(outfile, "%s\n", tokens[0].c_str());
      }
    }

    input.close();

    return Success;
}

}} // End namespaces

