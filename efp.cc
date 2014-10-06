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

#include "gamessparser.h"

INIT_PLUGIN

using namespace boost;
using namespace std;

namespace psi{ namespace psi4_gamess_efp {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "PSI4_GAMESS_MOS"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        options.add_str_i("GAMESS_OUTPUT_FILE", "");
    }

    return true;
}

extern "C" 
PsiReturnType psi4_gamess_efp(Options& options)
{
    int print = options.get_int("PRINT");

    boost::shared_ptr<PSIO> psio(_default_psio_lib_);
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");
    if(wfn->nirrep() != 1) throw PSIEXCEPTION("Must be run w/o symmetry for now.");
    boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));
    outfile->Printf("Nuclear repulsion energy = %20.12f\n", chkpt->rd_enuc());
    int nao, nso, nmo;
//    outfile->Printf("Number of AOs = %d\n", nao = chkpt->rd_nao());
//    outfile->Printf("Number of SOs = %d\n", nso = chkpt->rd_nso());
//    outfile->Printf("Number of MOs = %d\n", nmo = chkpt->rd_nmo());

//    SharedMatrix scf_SO = wfn->Ca_subset("SO", "ALL");
//    scf_SO->print(outfile);

    GamessOutputParser gparser(options);

    return Success;
}

}} // End namespaces

