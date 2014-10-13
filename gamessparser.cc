#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include "gamessparser.h"
#include "liboptions/liboptions.h"
#include "libmints/matrix.h"
#include "libmints/wavefunction.h"
#include "libmints/basisset.h"
#include "libpsio/psio.hpp"
#include "libchkpt/chkpt.hpp"
#include "psifiles.h"
#include <fstream>
#include <vector>


// A useful macro for defining spaces followed by a real number
#define SPACEFLOAT  "(?:\\s+(-?\\d+\\.\\d+))"

// The lines we're interested in look like this...
//  1  O  1  S    0.995825  -0.259033   0.005137   0.101188   0.000269
//  2  O  1  S    0.018925   1.035275  -0.038077  -0.767178  -0.004405
//  3  O  1 XX    0.000903   0.007263   0.001481   0.030810   0.000263
//  4  O  1 YY   -0.001200  -0.024227  -0.000911  -0.019347  -0.000197
//  5  O  1 ZZ    0.000297   0.016964  -0.000570  -0.011462  -0.000066
boost::regex coefs_re("^\\s+(\\d+)\\s+[A-Z]+\\s+\\d+\\s+[SXYZ]+" SPACEFLOAT SPACEFLOAT "?" SPACEFLOAT "?" SPACEFLOAT "?" SPACEFLOAT "?\\s*$");
boost::regex barenuc_re("^\\s+(\\d+)\\s+[A-Z]+\\s+\\d+\\s+[SXYZ]+" SPACEFLOAT SPACEFLOAT "?" SPACEFLOAT "?" SPACEFLOAT "?" SPACEFLOAT "?\\s*$");
boost::regex barenuc_end("^\\s+KINETIC ENERGY INTEGRALS\\s*$");

boost::regex nso_re("^\\s+NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =\\s*(\\d+)\\s*$");

// Regex objects to define regular expressions and capture match results
boost::smatch matchobj;

namespace psi{

void
GamessOutputParser::parse_H(std::ifstream &gamessout)
{
    if(nso_ == 0)
        throw PSIEXCEPTION("I should have the number of SOs by now...");

    HGamess_ = SharedMatrix(new Matrix("H from GAMESS", nso_, nso_));

    // zero indexed rows!
    int startrow = 0;
    int currow = 0;

    bool newblock = true;

    while(gamessout.good())
    {
        std::string line;
        std::getline(gamessout, line);
        // Look for something like the following...
        if (regex_match(line, matchobj, barenuc_re)){
            //int so = boost::lexical_cast<int>(matchobj[1])-1;
            newblock = false;
            for (int i = 2; i < matchobj.size(); ++i)
            {
                if(matchobj[i].length())
                {
                    double val;
                    try {val = boost::lexical_cast<double>(matchobj[i]);
                    }catch(...){
                        std::cout << "CANNOT CONVERT " << matchobj[i] << " TO DOUBLE";
                        exit(1);
                    }

                    (*HGamess_)(currow, startrow + i - 2) = val;
                    (*HGamess_)(startrow + i - 2, currow) = val;
                }
            }
            currow++;
        }
        else if(regex_match(line, barenuc_end))
            break;
        else if(!newblock)
        {
            // begin next block
            newblock = true;
            startrow += 5; // should be safe, even at the end, since we break right after
            currow = startrow;
        }
    }

#if 0
    std::vector< std::vector<double> >::const_iterator row_iter;
    std::vector<double>::const_iterator col_iter;
    outfile->Printf("\nH parsed from file, in GAMESS format\n");
    for(row_iter = gamess_H.begin(); row_iter < gamess_H.end(); ++row_iter){
        for(col_iter = row_iter->begin(); col_iter < row_iter->end(); ++col_iter){
            outfile->Printf("%9.6f ", *col_iter);
        }
        outfile->Printf("\n");
    }
#endif
}


void
GamessOutputParser::parse_mos(std::ifstream &gamessout)
{
    // Loop over the file, count MOs/SOs, and store coefficients
    std::vector< std::vector<double> > gamess_mos;
    int nmisses = 0;
    while(gamessout.good()){
        std::string line;
        std::getline(gamessout, line);
        // Look for something like the following...
        if (regex_match(line, matchobj, coefs_re)){
            nmisses = 0;
            std::vector<double> temp;
            int so = boost::lexical_cast<int>(matchobj[1])-1;
            for (int i = 2; i < matchobj.size(); ++i){
                if(matchobj[i].length()){
                    double val;
                    try {val = boost::lexical_cast<double>(matchobj[i]);
                    }catch(...){
                        std::cout << "CANNOT CONVERT " << matchobj[i] << " TO DOUBLE";
                        exit(1);
                    }
                    temp.push_back(val);
                }
            }
            if(so >= gamess_mos.size())
                // We haven't found any of these SOs yet; add this vector.
                gamess_mos.push_back(temp);
            else
                // We already have some MOs for this SO; append the current set.
                gamess_mos[so].insert(gamess_mos[so].end(), temp.begin(), temp.end());
        }else{
            // This isn't a matrix entry
            nmisses++;
        }
        if(nmisses > 5){
            // It's been 6 lines since we saw a valid entry; we're done.
            break;
        }
    }
#if 0
    outfile->Printf("\nMOs parsed from file, in GAMESS format\n");
    std::vector< std::vector<double> >::const_iterator row_iter;
    std::vector<double>::const_iterator col_iter;
    for(row_iter = gamess_mos.begin(); row_iter < gamess_mos.end(); ++row_iter){
        for(col_iter = row_iter->begin(); col_iter < row_iter->end(); ++col_iter){
            outfile->Printf("%9.6f ", *col_iter);
        }
        outfile->Printf("\n");
    }
#endif

    nmo_ = gamess_mos[0].size();

    if(gamess_mos.size() != nso_)
        throw PSIEXCEPTION("Error - Number of SOs doesn't match C!");

    CGamess_ = SharedMatrix(new Matrix("C from GAMESS", nso_, nmo_));
    for(int row = 0; row < nso_; ++row){
        for(int col = 0; col < nmo_; ++col){
            CGamess_->set(row, col, gamess_mos[row][col]);
        }
    }

}

void
GamessOutputParser::build_U_and_rotate()
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<BasisSet> basis = wfn->basisset();
    int nso = wfn->nso();
    int nmo = wfn->nmo();
    int nao = basis->nao();

    // A quick sanity check
    if(nao != nso_){
        std::cerr << "Expected " << nao << " AOs, but got " << nso_ << std::endl;
        exit(1);
    }
    if(nmo != nmo_){
        std::cerr << "Expected " << nmo << " MOs, but got " << nmo_ << std::endl;
        exit(1);
    }

    UC_ = SharedMatrix(new Matrix("C Transformation matrix", nso, nao));
    UH_ = SharedMatrix(new Matrix("H Transformation matrix", nso, nao));
    bool puream = basis->has_puream();
    int so_off = 0;
    int ao_off = 0;
    for(int shell = 0; shell < basis->nshell(); ++shell){
        int am = basis->shell(shell).am();
        /*
         * Re-map the SO indices.  The Cartesian ordering is as follows.
         *
         *         | 0 |   1   |         2         |                   3
         * --------|---|-------|-------------------|----------------------------------------
         *  Psi4:  | s | x y z | xx xy xz yy yz zz | xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
         *  GAMESS | s | x y z | xx yy zz xy xz yz | xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
         * --------|---|-------|-------------------|----------------------------------------
         *  offset | 0 | 0 1 2 |  0  3  4  1  5  2 |  0   3   4   5   9   7   1   6   8   2
         *
         * We also need to account for the missing normalization factor
         *      __________________________
         *     /         2^(2l)
         *    / --------------------------
         *  \/ (2lx-1)!!(2ly-1)!!(2lz-1)!!
         *
         * For sphericals, we need to do a little more.  It looks like GAMESS always reports
         * spherical harmonics in terms of their Cartesian components, so we just need to do
         * a quick remapping, as follows.  The Psi4 ordering is always 1, 1c, 1s, 2c, 2s,...
         *
         * l=1  | x | y | z
         * -----------------
         * 1_0  | 0 | 0 | 1 
         * 1_1c | 1 | 0 | 0 
         * 1_1s | 0 | 1 | 0
         *
         * l=2  | xx        | xy     | yy         | xz      | yz     | zz
         * --------------------------------------------------------------
         * 2_0  | -1/2      | 0      | -1/2       | 0       | 0      | 1 
         * 2_1c | 0         | 0      | 0          | sqrt(3) | 0      | 0 
         * 2_1s | 0         | 0      | 0          | 0       | sqrt(3)| 0 
         * 2_2c | sqrt(3)/2 | 0      | -sqrt(3)/2 | 0       | 0      | 0 
         * 2_2s | 0         | sqrt(3)| 0          | 0       | 0      | 0
         *
         * l=3 | xxx | xxy | xyy | yyy | xxz | xyz | yyz | xzz | yzz | zzz
         * ----------------------------------------------------------------
         * 3_0 | 0 | 0 | 0 | 0 | -frac32 | 0 | -frac32 | 0 | 0 | 1 
         * 3_1c| -sqrtfrac322 | 0 | -fracsqrtfrac322 | 0 | 0 | 0 | 0 | sqrt6 | 0 | 0 
         * 3_1s| 0 | -fracsqrtfrac322 | 0 | -fracsqrtfrac322 | 0 | 0 | 0 | 0 | sqrt6 | 0 
         * 3_2c| 0 | 0 | 0 | 0 | fracsqrt152 | 0 | -fracsqrt152 | 0 | 0 | 0 
         * 3_2s| 0 | 0 | 0 | 0 | 0 | sqrt15 | 0 | 0 | 0 | 0 
         * 3_3c| fracsqrtfrac522 | 0 | -frac3 sqrtfrac522 | 0 | 0 | 0 | 0 | 0 | 0 | 0 
         * 3_3s| 0 | frac3 sqrtfrac522 | 0 | -fracsqrtfrac522 | 0 | 0 | 0 | 0 | 0 | 0
         */
        if(puream){
            if(am == 0){
                // s
                UC_->set(so_off, ao_off, 1.0);
                UH_->set(so_off, ao_off, 1.0);
                so_off++; ao_off++;
            }else if(am == 1){
                // 1_0  <- z
                UC_->set(so_off+0, ao_off+2, 1.0);
                UH_->set(so_off+0, ao_off+2, 1.0);
                // 1_1c  <- x
                UC_->set(so_off+1, ao_off+0, 1.0);
                UH_->set(so_off+1, ao_off+0, 1.0);
                // 1_1s  <- y
                UC_->set(so_off+2, ao_off+1, 1.0);
                UH_->set(so_off+2, ao_off+1, 1.0);
                so_off += 3; ao_off += 3;
            }else if(am == 2){
                double Caa = 1.0;
                double Cab = sqrt(3.0);
                double Haa = 1.0;
                double Hab = 1.0/sqrt(3.0);

                double C0 = 2.0/3.0;
                double C1c = 1.0/sqrt(3.0);
                double C1s = 1.0/sqrt(3.0);
                double C2c = 1.0/sqrt(3.0);
                double C2s = 1.0/sqrt(3.0);
                double H0 = 1.0;
                double H1c = sqrt(3.0);
                double H1s = sqrt(3.0);
                double H2c = sqrt(3.0)/2.0;
                double H2s = sqrt(3.0);
                // 2_0 <- zz - 0.5xx - 0.5yy
                UC_->set(so_off+0, ao_off+2, 1.0*C0*Caa);
                UC_->set(so_off+0, ao_off+0, -0.5*C0*Caa);
                UC_->set(so_off+0, ao_off+1, -0.5*C0*Caa);
                UH_->set(so_off+0, ao_off+2, 1.0*H0*Haa);
                UH_->set(so_off+0, ao_off+0, -0.5*H0*Haa);
                UH_->set(so_off+0, ao_off+1, -0.5*H0*Haa);
                // 2_1c <- xz
                UC_->set(so_off+1, ao_off+4, 1.0*C1c*Cab);
                UH_->set(so_off+1, ao_off+4, 1.0*H1c*Hab);
                // 2_1s <- yz
                UC_->set(so_off+2, ao_off+5, 1.0*C1s*Cab);
                UH_->set(so_off+2, ao_off+5, 1.0*H1s*Hab);
                // 2_2c <- xx - yy
                UC_->set(so_off+3, ao_off+0, 1.0*C2c*Haa);
                UC_->set(so_off+3, ao_off+1, -1.0*C2c*Haa);
                UH_->set(so_off+3, ao_off+0, 1.0*H2c*Haa);
                UH_->set(so_off+3, ao_off+1, -1.0*H2c*Haa);
                // 2_2s <- xy
                UC_->set(so_off+4, ao_off+3, 1.0*C2s*Cab);
                UH_->set(so_off+4, ao_off+3, 1.0*H2s*Hab);
                so_off += 5; ao_off += 6;
            }else if(am == 3){
                double Caaa = 1.0;
                double Caab = sqrt(5.0);
                double Cabc = sqrt(15.0);
                double Haaa = 1.0;
                double Haab = 1.0/sqrt(5.0);
                double Habc = 1.0/sqrt(15.0);

                double C0 = 2.0/11.0; // good
                double C1c = 1.0/sqrt(15.0); // bad
                double C1s = 1.0/sqrt(15.0); // bad
                double C2c = 1.0/sqrt(15.0); // good
                double C2s = 1.0/sqrt(15.0); // good
                double C3c = 1.0/sqrt(15.0); // bad
                double C3s = 1.0/sqrt(15.0); // bad
                double H0 = 1.0;
                double H1c = sqrt(6.0);
                double H1s = sqrt(6.0);
                double H2c = sqrt(15.0/4.0);
                double H2s = sqrt(15.0);
                double H3c = sqrt(5.0/2.0);
                double H3s = sqrt(5.0/2.0);
         //   xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
         //    0   1   2   3   4   5   6   7   8   9
                // 3_0 <- zzz - 1.5xxz - 1.5yyz
                UC_->set(so_off+0, ao_off+2, 1.0*C0*Caaa);
                UC_->set(so_off+0, ao_off+4, -1.5*C0*Caab);
                UC_->set(so_off+0, ao_off+6, -1.5*C0*Caab);
                UH_->set(so_off+0, ao_off+2, 1.0*H0*Haaa);
                UH_->set(so_off+0, ao_off+4, -1.5*H0*Haab);
                UH_->set(so_off+0, ao_off+6, -1.5*H0*Haab);
                // 3_1c <- sqrt(6)xzz - sqrt(3/8)xyy - sqrt(3/8)xxx
                UC_->set(so_off+1, ao_off+7, 1.0*C1c*Caab);
                UC_->set(so_off+1, ao_off+5, -0.25*C1c*Caab);
                UC_->set(so_off+1, ao_off+0, -0.25*C1c*Caaa);
                UH_->set(so_off+1, ao_off+7, 1.0*H1c*Haab);
                UH_->set(so_off+1, ao_off+5, -0.25*H1c*Haab);
                UH_->set(so_off+1, ao_off+0, -0.25*H1c*Haaa);
                // 3_1s <- sqrt(6)yzz - sqrt(3/8)yyy - sqrt(3/8)xxy
                UC_->set(so_off+2, ao_off+8, 1.0*C1s*Caab);
                UC_->set(so_off+2, ao_off+1, -0.25*C1s*Caaa);
                UC_->set(so_off+2, ao_off+3, -0.25*C1s*Caab);
                UH_->set(so_off+2, ao_off+8, 1.0*H1s*Haab);
                UH_->set(so_off+2, ao_off+1, -0.25*H1s*Haaa);
                UH_->set(so_off+2, ao_off+3, -0.25*H1s*Haab);
                // 3_2c <- sqrt(15/4)xxz - sqrt(15/4)yyz
                UC_->set(so_off+3, ao_off+4, 1.0*C2c*Caab);
                UC_->set(so_off+3, ao_off+6, -1.0*C2c*Caab);
                UH_->set(so_off+3, ao_off+4, 1.0*H2c*Haab);
                UH_->set(so_off+3, ao_off+6, -1.0*H2c*Haab);
                // 3_2s <- sqrt(15)xyz
                UC_->set(so_off+4, ao_off+9, 1.0*C2s*Cabc);
                UH_->set(so_off+4, ao_off+9, 1.0*H2s*Habc);
                // 3_3c <- sqrt(5/2)/2xxx - 3sqrt(5/2)/2xyy
                UC_->set(so_off+5, ao_off+0, 0.5*C3c*Caaa);
                UC_->set(so_off+5, ao_off+5, -1.5*C3c*Caab);
                UH_->set(so_off+5, ao_off+0, 0.5*H3c*Haaa);
                UH_->set(so_off+5, ao_off+5, -1.5*H3c*Haab);
                // 3_3s <- 3sqrt(5/2)/2xxy - sqrt(5/2)/2yyy
                UC_->set(so_off+6, ao_off+3, 1.5*C3s*Caab);
                UC_->set(so_off+6, ao_off+1, -0.5*C3s*Caaa);
                UH_->set(so_off+6, ao_off+3, 1.5*H3s*Haab);
                UH_->set(so_off+6, ao_off+1, -0.5*H3s*Haaa);
                so_off += 7; ao_off += 10;
            }else{
                throw PSIEXCEPTION("f functions not yet implemented for pure A.M.");
            }
        }else{
            if(am == 0){
                // s
                UC_->set(so_off, ao_off, 1.0);
                UH_->set(so_off, ao_off, 1.0);
                so_off++; ao_off++;
            }else if(am == 1){
                // x
                UC_->set(so_off+0, ao_off+0, 1.0);
                UH_->set(so_off+0, ao_off+0, 1.0);
                // y
                UC_->set(so_off+1, ao_off+1, 1.0);
                UH_->set(so_off+1, ao_off+1, 1.0);
                // z
                UC_->set(so_off+2, ao_off+2, 1.0);
                UH_->set(so_off+2, ao_off+2, 1.0);
                so_off += 3; ao_off += 3;
            }else if(am == 2){
                double Caa = 1.0;
                double Cab = sqrt(3.0);
                double Haa = 1.0;
                double Hab = 1.0/sqrt(3.0);
                // xx
                UC_->set(so_off+0, ao_off+0, Caa);
                UH_->set(so_off+0, ao_off+0, Haa);
                // xy
                UC_->set(so_off+1, ao_off+3, Cab);
                UH_->set(so_off+1, ao_off+3, Hab);
                // xz
                UC_->set(so_off+2, ao_off+4, Cab);
                UH_->set(so_off+2, ao_off+4, Hab);
                // yy
                UC_->set(so_off+3, ao_off+1, Caa);
                UH_->set(so_off+3, ao_off+1, Haa);
                // yz
                UC_->set(so_off+4, ao_off+5, Cab);
                UH_->set(so_off+4, ao_off+5, Hab);
                // zz
                UC_->set(so_off+5, ao_off+2, Caa);
                UH_->set(so_off+5, ao_off+2, Haa);
                so_off += 6; ao_off += 6;
            }else if(am == 3){
                double Caaa = 1.0;
                double Caab = sqrt(5.0);
                double Cabc = sqrt(15.0);
                double Haaa = 1.0;
                double Haab = 1.0/sqrt(5.0);
                double Habc = 1.0/sqrt(15.0);
                // xxx
                UC_->set(so_off+0, ao_off+0, Caaa);
                UH_->set(so_off+0, ao_off+0, Haaa);
                // xxy
                UC_->set(so_off+1, ao_off+3, Caab);
                UH_->set(so_off+1, ao_off+3, Haab);
                // xxz
                UC_->set(so_off+2, ao_off+4, Caab);
                UH_->set(so_off+2, ao_off+4, Haab);
                // xyy
                UC_->set(so_off+3, ao_off+5, Caab);
                UH_->set(so_off+3, ao_off+5, Haab);
                // xyz
                UC_->set(so_off+4, ao_off+9, Cabc);
                UH_->set(so_off+4, ao_off+9, Habc);
                // xzz
                UC_->set(so_off+5, ao_off+7, Caab);
                UH_->set(so_off+5, ao_off+7, Haab);
                // yyy
                UC_->set(so_off+6, ao_off+1, Caaa);
                UH_->set(so_off+6, ao_off+1, Haaa);
                // yyz
                UC_->set(so_off+7, ao_off+6, Caab);
                UH_->set(so_off+7, ao_off+6, Haab);
                // yzz
                UC_->set(so_off+8, ao_off+8, Caab);
                UH_->set(so_off+8, ao_off+8, Haab);
                // zzz
                UC_->set(so_off+9, ao_off+2, Caaa);
                UH_->set(so_off+9, ao_off+2, Haaa);
                so_off += 10; ao_off += 10;
            }else{
                throw PSIEXCEPTION("f functions not yet implemented for pure A.M.");
            }
        }
    }

    C_ = SharedMatrix(new Matrix("C reordered for Psi4", nso, nmo));
    // C = U Cgamess
    C_->gemm(false, false, 1.0, UC_, CGamess_, 0.0);

    if(H_found_)
    {
        H_ = SharedMatrix(new Matrix("H reordered for Psi4", nso, nso));
        // H = Ut Hgamess U
        H_->back_transform(HGamess_, UH_);
    }

    if(options_.get_int("PRINT") > 1){
        wfn->Ca()->print();
        CGamess_->print();
        UC_->print();
        C_->print();

        if(H_found_)
        {
            HGamess_->print();
            UH_->print();
            H_->print();
            wfn->H()->print();
        }


    }

    wfn->Ca()->copy(C_);
    wfn->Cb()->copy(C_);

    wfn->save();

    //exit(1);
}


GamessOutputParser::GamessOutputParser(Options &options):
    options_(options), H_found_(false)
{
    nso_ = 0;

    // The gamess output file
    std::ifstream gamessout(options.get_str("GAMESS_OUTPUT_FILE").c_str());
    if (!gamessout)
        throw PSIEXCEPTION("Unable to open the GAMESS output file.");

    bool mos_found = false;

    // A line with only spaces and "EIGENVECTORS" marks the start of the MO coefficients
    boost::regex evecs_label_re("^\\s+EIGENVECTORS\\s*$");
    boost::regex H_label_re("BARE NUCLEUS HAMILTONIAN INTEGRALS");
    while(gamessout.good()){
        std::string line;
        std::getline(gamessout, line);
        if(regex_match(line, matchobj, nso_re))
        {
            try
            {
                nso_ = boost::lexical_cast<int>(matchobj[1]);
                std::cout << "====NSO: " << nso_ << "\n";
            }
            catch(...)
            {
                        std::cout << "CANNOT CONVERT " << matchobj[1] << " TO INT";
                        exit(1);
            }
        }
        // Look for the start of the H matrix definition
        if (regex_match(line, matchobj, evecs_label_re)){
            mos_found = true;
            parse_mos(gamessout);
        }
        // Look for the start of the H definition
        if (regex_search(line, matchobj, H_label_re)){
            H_found_ = true;
            parse_H(gamessout);
        }
    }

    if(!H_found_)
        outfile->Printf("WARNING - No H matrix was found in the GAMESS output file provided\n");

    if(!mos_found)
        throw PSIEXCEPTION("No MOs were found in the GAMESS output file provided");
    build_U_and_rotate();
}

} // End namespace
