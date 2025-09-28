#ifndef TREEMANAGER_h
#define TREEMANAGER_h

#include "FaserCalRecoPlot.h"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include <string>
class TreeManager {
        public :
                TreeManager(std::string infilename);
                virtual ~TreeManager();
                
                FaserCalRecoPlot *tmCD()     { return m_tres; };
        protected:
                void          init();
                void          nullify();
                
        private:
                std::string   m_inFilename;
                
                TFile         *m_inFile;
                
                TTree         *m_treeRes;
                
                FaserCalRecoPlot *m_tres;
};
#endif //TreeManager_h
