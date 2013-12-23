#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libciomr/libciomr.h>
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

    double **scf_old = chkpt->rd_scf();
//    fprintf(outfile, "SCF Matrix from Checkpoint File:\n");
//    print_mat(scf_old, nso, nmo, outfile);
    free_block(scf_old);

    SharedMatrix scf_SO = wfn->Ca_subset("SO", "ALL"); // AO vs. SO definition in PSI4 is screwed up
    scf_SO->print(outfile);

    boost::shared_ptr<BasisSet> basisset = wfn->basisset();
    boost::shared_ptr<IntegralFactory> integral = wfn->integral();
    PetiteList petite(basisset, integral, true);
    SharedMatrix aotoso = petite.aotoso();
    // aotoso->print(outfile);

    SharedMatrix scf_AO(new Matrix("SCF AO x MO", nao, nmo));
    scf_AO->gemm(false, false, 1.0, aotoso, scf_SO, 0.0);
    scf_AO->print(outfile);   

    // Read MO data from GAMESS output file
    ifstream input;
    input.open("mos.txt");
    if(!input.good()) throw PSIEXCEPTION("Error opening mos.txt file.");
    char buf[512];
    SharedVector energies = SharedVector(new Vector("MO Energies", nmo));
    int line_count = 0;
    int block_count = 0;
    while(!input.eof()) {
      input.getline(buf, 512);  // grab a line from input
      stringstream cppbuf(buf);
      string token;
      vector<string> tokens;
      while(cppbuf >> token) tokens.push_back(token);
      int mo_min = block_count * 5;
      int mo_max = mo_min + ((block_count == nmo/5) ? (nmo % mo_min) : 4);

      if(line_count == 1) { // Grab MO energies
        for(int i=mo_min; i <= mo_max; i++)
          energies->set(i,atof(tokens[i % 5].c_str()));
        line_count++;
      }
      else if(line_count == (nao+3)) { block_count++; line_count++; }
      else if(line_count > 2 && line_count < (3+nao)) {
        for(int i=mo_min; i <= mo_max; i++)
          scf_AO->set(line_count-3, i, atof(tokens[(i % 5) + 4].c_str()));
        line_count++;
      }
      else line_count++;
      // skip blank lines and mark start of data-block
      if(!tokens.size()) line_count = 0;
    }
 
    energies->print(outfile);
    chkpt->wt_evals(energies->pointer());
//    scf_AO->print(outfile);

    input.close();

    // Transform scf_AO to scf_SO
//    SharedMatrix scf_SO(new Matrix("SCF SO x MO", nso, nmo));
//    scf_SO->gemm(true, false, 1.0, aotoso, scf_AO, 0.0);
//    scf_SO->print(outfile);

    double **p_scf_SO = scf_SO->pointer();
    chkpt->wt_scf(p_scf_SO);   
    chkpt->wt_alpha_scf(p_scf_SO);   

    // Check to make sure it worked
    double **scf_new = chkpt->rd_alpha_scf();
//    fprintf(outfile, "SCF Matrix from Checkpoint File:\n");
//    print_mat(scf_new, nso, nmo, outfile);
    free_block(scf_new);

    return Success;
}

}} // End namespaces

