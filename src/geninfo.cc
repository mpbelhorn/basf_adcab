// ******************************************************************
// Generator Information functions.
// by A. Zupanc
//
// ******************************************************************

#include "belle.h"
#include "panther/panther.h"
#include "mdst/mdst.h"
#include MDST_H
#include HEPEVT_H

#include "particle/utility.h"
#include "geninfo.h"
#include "userinfo.h"

//#include "user_def.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;

int isFSP(Gen_hepevt P)
{
  switch (abs(P.idhep())) {
    case 211:
      return 1;
    case 321:
      return 1;
    case 11:
      return 1;
    case 13:
      return 1;
    case 22:
      return 1;
    case 2212:
      return 1;
    case 111:
      return 1;
    case 310:
      return 1;
    case 130:
      return 1;
    case 2112:
      return 1;
    default:
      return 0;
  }
}

void appendRecFSP(Particle p, std::vector<Particle> &children)
{
  // std::cout << "Particle Decay mode = "
  //           << dynamic_cast<UserInfo&>(p.userInfo()).decayMode()
  //           << std::endl;
  for (int i = 0; i < (int)p.nChildren(); ++i) {
    if (p.child(i).nChildren() &&
        p.child(i).lund() != 111 &&
        p.child(i).lund()!=310) {
      appendRecFSP(p.child(i),children);
    } else {
      // std::cout << "Appending FSP: " << p.child(i).lund() << std::endl;
      children.push_back(p.child(i));
    }
  }
}

int appendGenFSP(const Gen_hepevt &gen, std::vector<Gen_hepevt> &children)
{
  Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();

  // int ndaug = (gen.daLast()-gen.daFirst()) + 1;
  // std::cout << "Appending FS children for "
  //           << gen.idhep() << " with nChildren = " << ndaug << std::endl;
  for (int i = gen.daFirst(); i <= gen.daLast(); ++i) {
    //std::cout << "searching for particle with ID = " << i << std::endl;
    if (i == 0) {
      std::cout << "[Zupanc] appendGenFSP: "
                << "requesting Particle with ID = 0! Exiting."
                << std::endl;
      return 0;
    }
    const Gen_hepevt& child = GenMgr(Panther_ID(i));
    // if (!child) {
    //   std::cout << "Particle with this PantherID does not exist."
    //             << std::endl;
    // }
    int ndaug2 = (child.daLast()-child.daFirst()) + 1;
    // std::cout << " -> child " << i << ": ID = " << i
    //           << "with nChildren = " << ndaug2 << std::endl;
    if (ndaug2 && !isFSP(child)) {
      // std::cout << "not FSP: " << child.idhep() << std::endl;
      appendGenFSP(child, children);
    } else {
      // std::cout << "Appending FSP: " << child.idhep() << std::endl;
      children.push_back(child);
    }
  }
  return 1;
}

int dumpGenFSP(const Gen_hepevt &gen)
{
  Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();

  int ndaug = (gen.daLast()-gen.daFirst()) + 1;
  std::cout << "Appending FS children for " << gen.idhep()
            << " with nChildren = " << ndaug << std::endl;
  for (int i = gen.daFirst(); i <= gen.daLast(); ++i) {
    std::cout << "searching for particle with ID = " << i << std::endl;
    if (i == 0) {
      std::cout << "[Zupanc] appendGenFSP: "
                << "requesting Particle with ID = 0! Exiting." << std::endl;
      return 0;
    }
    const Gen_hepevt& child = GenMgr(Panther_ID(i));
    if (!child) {
      std::cout << "Particle with this PantherID does not exist." << std::endl;
    }
    int ndaug2 = (child.daLast()-child.daFirst()) + 1;
    std::cout << " -> child " << i << ": ID = " << i
              << "with nChildren = " << ndaug2 << std::endl;
    if (ndaug2 && !isFSP(child)) {
      std::cout << "not FSP: " << child.idhep() << std::endl;
      dumpGenFSP(child);
    } else {
      std::cout << "Appending FSP: " << child.idhep() << std::endl;
    }
  }
  return 1;
}

int commonMother(std::vector<int> &mothers)
{
  // std::cout << "...common mother of all..." << std::endl;
  if (mothers.size() == 0) {
    // std::cout << "size of mothers vector == 0! Exit." << std::endl;
    return 0;
  } else if (mothers.size() == 1) {
    // std::cout << "size of mothers vector == 1! Return mother = "
    //           << mothers[0] << std::endl;
    return mothers[0];
  }
  int motherID = mothers[0];
  // std::cout << "Searching for ID = " << motherID << std::endl;
  for (int i = 1; i < (int)mothers.size(); ++i) {
    //std::cout << "[" << i << "] ID = " << mothers[i] << std::endl;
    if (motherID!=mothers[i]) {
      //std::cout << "Different IDs! Exit." << std::endl;
      return 0;
    }
  }
  return motherID;
}

void fillMothers(Particle &p, std::vector<int> &thisMothers,
    std::vector<int> &otherMothers)
{
  Gen_hepevt motherThis(p.child(0).genHepevt());

  // std::cout << "This particle = " << p.child(0).genHepevt().idhep()
  //           << " with ID = " << p.child(0).genHepevt().get_ID()
  //           << std::endl;
  while(motherThis.mother()) {
    motherThis = motherThis.mother();
    // std::cout << " -> mother = " << motherThis.idhep()
    //           << " ID = " << motherThis.get_ID() << std::endl;
    thisMothers.push_back(motherThis.get_ID());
  }

  for(int i = 1; i<(int)p.nChildren(); ++i) {
    motherThis = p.child(i).genHepevt();
    // std::cout << "Other particle (" << i << ") = "
    //           << motherThis.idhep() << " with ID = "
    //           << motherThis.get_ID() << std::endl;
    while(motherThis.mother()) {
      motherThis = motherThis.mother();
      // std::cout << " -> mother = " << motherThis.idhep()
      //           << " ID = " << motherThis.get_ID() << std::endl;
      otherMothers.push_back(motherThis.get_ID());
    }
  }
}

void fillMothers(Particle &A, Particle &B,
    std::vector<int> &thisMothers, std::vector<int> &otherMothers)
{
  Gen_hepevt motherThis(A.genHepevt());
  while(motherThis.mother()) {
    motherThis = motherThis.mother();
    thisMothers.push_back(motherThis.get_ID());
  }

  motherThis = B.genHepevt();
  while(motherThis.mother()) {
    motherThis = motherThis.mother();
    otherMothers.push_back(motherThis.get_ID());
  }
}

void fillMothers(Particle &A, Particle &B, Particle &C,
    std::vector<int> &thisMothers, std::vector<int> &otherMothers)
{
  Gen_hepevt motherThis(A.genHepevt());

  while(motherThis.mother()) {
    motherThis = motherThis.mother();
    thisMothers.push_back(motherThis.get_ID());
  }

  motherThis = B.genHepevt();
  while(motherThis.mother()) {
    motherThis = motherThis.mother();
    otherMothers.push_back(motherThis.get_ID());
  }

  motherThis = C.genHepevt();
  while(motherThis.mother()) {
    motherThis = motherThis.mother();
    otherMothers.push_back(motherThis.get_ID());
  }
}

int getCommonMother(Particle &A, Particle &B)
{
  if (! &A.userInfo() ) setUserInfo(A);
  if (! &B.userInfo() ) setUserInfo(B);

  if (!(dynamic_cast<UserInfo&>(A.userInfo()).genHepevtChecked())) {
    setMCtruth(A);
  }
  if (!(dynamic_cast<UserInfo&>(B.userInfo()).genHepevtChecked())) {
    setMCtruth(B);
  }

  if (!A.genHepevt())
    return 0;
  if (!B.genHepevt())
    return 0;

  std::vector<int> thisMothers;
  std::vector<int> otherMothers;

  fillMothers(A, B, thisMothers, otherMothers);
  if (thisMothers.size() == 0 || otherMothers.size() == 0) {
    return 0;
  }

  int motherID = findCommonMother(2, thisMothers, otherMothers);

  if (!motherID) {
    std::cout << "[Zupanc] setMCtruth: "
              << "requesting particle with ID = 0! Exiting." << std::endl;
    return 0;
  }

  Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
  const Gen_hepevt& thisGen = GenMgr(Panther_ID(motherID));

  return thisGen.idhep();
}

int getCommonMother(Particle &A, Particle &B, Particle &C) {
  if (! &A.userInfo()) setUserInfo(A);
  if (! &B.userInfo()) setUserInfo(B);
  if (! &C.userInfo()) setUserInfo(C);

  if (!(dynamic_cast<UserInfo&>(A.userInfo()).genHepevtChecked())) {
    setMCtruth(A);
  }
  if (!(dynamic_cast<UserInfo&>(B.userInfo()).genHepevtChecked())) {
    setMCtruth(B);
  }
  if (!(dynamic_cast<UserInfo&>(C.userInfo()).genHepevtChecked())) {
    setMCtruth(C);
  }
  if (!A.genHepevt()) return 0;
  if (!B.genHepevt()) return 0;
  if (!C.genHepevt()) return 0;

  std::vector<int> thisMothers;
  std::vector<int> otherMothers;

  fillMothers(A, B, C, thisMothers, otherMothers);
  if (thisMothers.size() == 0 || otherMothers.size() == 0) {
    return 0;
  }

  int motherID = findCommonMother(3, thisMothers, otherMothers);

  if (!motherID) {
    std::cout << "[Zupanc] setMCtruth: "
              << "requesting particle with ID = 0! Exiting." << std::endl;
    return 0;
  }

  Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
  const Gen_hepevt& thisGen = GenMgr(Panther_ID(motherID));

  return thisGen.idhep();
}

int findCommonMother(int nChildren,
    std::vector<int> thisMothers, std::vector<int> otherMothers) {
  for (int i = 0; i < (int)thisMothers.size(); ++i) {
    int counter = 0;
    for (int j = 0; j < (int)otherMothers.size(); ++j) {
      if (thisMothers[i] == otherMothers[j]) counter++;
    }
    if (counter == nChildren-1) return thisMothers[i];
  }
  return 0;
}

int findCommonMother(Gen_hepevt pThis, Gen_hepevt pOther, int level) {
  Gen_hepevt motherThis(pThis);

  // std::cout << "This = " << pThis.idhep() << " [ID = " << pThis.get_ID()
  //           << "] vs. " << "Other = " << pOther.idhep() << " [ID = "
  //           << pOther.get_ID() << std::endl;
  int i = 1;
  while(motherThis.mother()) {
    motherThis = motherThis.mother();
    // std::cout << "[" << i << "] This mother = " << motherThis.idhep()
    //           << " ID = " << motherThis.get_ID() << std::endl;
    i++;

    if (i > level) {
      // std::cout << "Other: " << pOther.idhep()
      //           << " [ID = " << pOther.get_ID() << std::endl;

      Gen_hepevt motherOther(pOther);
      while(motherOther.mother()) {
        motherOther = motherOther.mother();
        // std::cout << "[" << i << "] Other mother = " << motherOther.idhep()
        //           << " ID = " << motherOther.get_ID() << std::endl;
        if (motherThis.get_ID() == motherOther.get_ID()) {
          return motherThis.get_ID();
        }
      }
    }
  }
  return 0;
}

int compareFinalStates(std::vector<Particle> reconstructed,
    std::vector<Gen_hepevt> generated)
{
  if (reconstructed.size() == generated.size()) {
    int missID = 0;
    int missPi0 = 0;
    for (int i = 0; i < (int)reconstructed.size(); ++i) {
      if (reconstructed[i].genHepevt()) {
        int link = 0;
        for (int j = 0; j < (int)generated.size(); ++j) {
          if (reconstructed[i].genHepevt().get_ID() == generated[j].get_ID()) {
            link = 1;
            if (reconstructed[i].lund()!=generated[j].idhep())
              missID++;
            if (reconstructed[i].lund() == 111 &&
                dynamic_cast<UserInfo&>(reconstructed[i].userInfo()).genHepevtLink() == 5) {
              missPi0++;
            }
            break;
          }
        }
        if (!link) {
          std::cout << "[Zupanc] compareFinalStates: "
                    << "Particle (" << reconstructed[i].lund()
                    << " with hepevt = " << reconstructed[i].genHepevt().idhep()
                    << ") not found in list of FSP (gen)!" << std::endl;
          return -11;
        }
      } else {
        std::cout << "[Zupanc] compareFinalStates: "
                  << "Particle without link to genHepevt! [-10]" << std::endl;
        return -10;
      }
    }
    if (missID && missPi0) return 6;
    else if (missID) return 2;
    else if (missPi0) return 5;
    return 1;
  } else if (reconstructed.size() < generated.size()) { // missing particle
    int missing = 0;
    std::vector<int> missP;
    for (int i = 0; i < (int)generated.size(); ++i) {
      int link = 0;
      for (int j = 0; j < (int)reconstructed.size(); ++j) {
        if (reconstructed[j].genHepevt()) {
          if (reconstructed[j].genHepevt().get_ID() == generated[i].get_ID()) {
            link = 1;
            break;
          }
        } else {
          std::cout << "[Zupanc] compareFinalStates: "
                    << "Particle without link to genHepevt! [-9]" << std::endl;
          return -9;
        }
      }
      if (!link) {
        missing++;
        missP.push_back(i);
      }
    }
    if (missing) {
      int masslessOnly = 1;
      for (int i = 0; i < (int)missP.size(); ++i) {
        if (generated[missP[i]].idhep() != 22 &&
            abs(generated[missP[i]].idhep()) != 12 &&
            abs(generated[missP[i]].idhep())!=14 &&
            abs(generated[missP[i]].idhep())!=16) {
          masslessOnly = 0;
        }
      }
      if (masslessOnly) {
        return 4;
      } else {
        return 3;
      }
    } else {
      std::cout << "[Zupanc] compareFinalStates: "
                << "At least one particle should be missing!" << std::endl;
      return -8;
    }
  } else {
    std::cout << "[Zupanc] compareFinalStates: "
              << "More reconstructed than generated particle in final state!"
              << std::endl;
    return -5;
  }
}

void setMCtruthPi0(Particle &p)
{
  if (! &p.userInfo()) setUserInfo(p); // create UserInfo, if it does not exsist

  if (dynamic_cast<UserInfo&>(p.userInfo()).genHepevtChecked() == 1) {
    //std::cout << "Nothing to do - particle already checked!" << std::endl;
    return;
  }

  // set checked flag to 1
  dynamic_cast<UserInfo&>(p.userInfo()).genHepevtChecked(1);
  //std::cout << "[Zupanc] setMCtruthPi0: start" << std::endl;
  for(int i=0; i<(int)p.nChildren(); ++i) {
    if (! &p.child(i).userInfo()) setUserInfo(p.child(i));
    if (p.child(i).mdstGamma()) {
      const Gen_hepevt & hep(get_hepevt(p.child(i).mdstGamma()));
      dynamic_cast<UserInfo&>(p.child(i).userInfo()).genHepevtChecked(1);
      if (hep) {
        p.child(i).relation().genHepevt(hep);
        dynamic_cast<UserInfo&>(p.child(i).userInfo()).genHepevtLink(1);
        // std::cout << " link found ("
        //           << dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink()
        //           << ") " << std::endl;
      } else {
        std::cout << " [Zupanc] setMCtruthPi0: child "
                  << i << " has no link! " << std::endl;
      }
    } else {
      std::cout << " [Zupanc] setMCtruthPi0: pi0 child is not mdstGamma! Exit."
                << std::endl;
    }
  }

  Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
  if (p.child(0).genHepevt() && p.child(1).genHepevt()) {
    const Gen_hepevt & mother1(p.child(0).genHepevt());
    const Gen_hepevt & mother2(p.child(1).genHepevt());
    if (mother1.mother() && mother2.mother()) {
      if (mother1.mother().get_ID() == mother2.mother().get_ID()) {
        if (mother1.mother().get_ID() == 0) {
          std::cout << "[Zupanc] setMCtruthPi0: "
                    << "requesting particle with ID = 0! Exiting." << std::endl;
          return;
        }
        const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother1.mother().get_ID()));
        p.relation().genHepevt(thisGen);
        dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(1);
        return;
      } else if (mother1.mother().idhep() == 111 &&
          mother2.mother().idhep() == 111) {
        // std::cout << "[Zupanc] setMCtruthPi0: pi0 different mothers "
        //           << mother1.mother().idhep() << " vs "
        //           << mother2.mother().idhep() << std::endl;
        if (p.child(0).e() > p.child(1).e()) {
          if (mother1.mother().get_ID() == 0) {
            std::cout << "[Zupanc] setMCtruthPi0: "
                      << "requesting particle with ID = 0! Exiting."
                      << std::endl;
            return;
          }
          const Gen_hepevt& thisGen = GenMgr(
              Panther_ID(mother1.mother().get_ID()));
          p.relation().genHepevt(thisGen);
          dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(4);
          return;
        } else {
          if (mother2.mother().get_ID() == 0) {
            std::cout << "[Zupanc] setMCtruthPi0: "
                      << "requesting particle with ID = 0! Exiting."
                      << std::endl;
            return;
          }
          const Gen_hepevt& thisGen = GenMgr(
              Panther_ID(mother2.mother().get_ID()));
          p.relation().genHepevt(thisGen);
          dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(4);
          return;
        }
      } else if (mother1.mother().idhep() == 111 ||
          mother2.mother().idhep() == 111) {
        if (mother1.mother().idhep() == 111) {
          if (mother1.mother().get_ID() == 0) {
            std::cout << "[Zupanc] setMCtruthPi0: "
                      << "requesting particle with ID = 0! Exiting."
                      << std::endl;
            return;
          }
          const Gen_hepevt& thisGen = GenMgr(
              Panther_ID(mother1.mother().get_ID()));
          p.relation().genHepevt(thisGen);
          dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(3);
          return;
        } else {
          if (mother2.mother().get_ID() == 0) {
            std::cout << "[Zupanc] setMCtruthPi0: "
                      << "requesting particle with ID = 0! Exiting."
                      << std::endl;
            return;
          }
          const Gen_hepevt& thisGen = GenMgr(
              Panther_ID(mother2.mother().get_ID()));
          p.relation().genHepevt(thisGen);
          dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(3);
          return;
        }
      }
    } else {
      // std::cout << "[Zupanc] setMCtruthPi0: "
      //           << "one child without mother " << std::endl;
      if (mother1.mother()) {
        if (mother1.mother().get_ID() == 0) {
          std::cout << "[Zupanc] setMCtruthPi0: "
                    << "requesting particle with ID = 0! Exiting."
                    << std::endl;
          return;
        }
        const Gen_hepevt& thisGen = GenMgr(
            Panther_ID(mother1.mother().get_ID()));
        p.relation().genHepevt(thisGen);
        dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(2);
        return;
      } else if (mother2.mother()) {
        if (mother2.mother().get_ID() == 0) {
          std::cout << "[Zupanc] setMCtruthPi0: "
                    << "requesting particle with ID = 0! Exiting."
                    << std::endl;
          return;
        }
        const Gen_hepevt& thisGen = GenMgr(
            Panther_ID(mother2.mother().get_ID()));
        p.relation().genHepevt(thisGen);
        dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(2);
        return;
      }
    }
  } else if (p.child(0).genHepevt()) {
    std::cout << "[Zupanc] setMCtruthPi0: "
              << "Only first photon with link to gen_hepevt: " << std::endl;
    const Gen_hepevt & mother(p.child(0).genHepevt());
    std::cout << " -> idhep = " << mother.idhep() << std::endl;
    if (mother.mother()) {
      std::cout << " -> mother ID = " << mother.mother().idhep() << std::endl;
      if (mother.mother().idhep() == 111) {
        if (mother.mother().get_ID() == 0) {
          std::cout << "[Zupanc] setMCtruthPi0: "
                    << "requesting particle with ID = 0! Exiting."
                    << std::endl;
          return;
        }
        const Gen_hepevt& thisGen = GenMgr(
            Panther_ID(mother.mother().get_ID()));
        p.relation().genHepevt(thisGen);
        dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(5);
        return;
      }
      return;
    }
  } else if (p.child(1).genHepevt()) {
    std::cout << "[Zupanc] setMCtruthPi0: "
              << "Only second photon with link to gen_hepevt: "
              << std::endl;
    const Gen_hepevt & mother(p.child(1).genHepevt());
    std::cout << " -> idhep = " << mother.idhep() << std::endl;
    if (mother.mother()) {
      std::cout << " -> mother ID = " << mother.mother().idhep() << std::endl;
      if (mother.mother().idhep() == 111) {
        if (mother.mother().get_ID() == 0) {
          std::cout << "[Zupanc] setMCtruthPi0: "
                    << "requesting particle with ID = 0! Exiting."
                    << std::endl;
          return;
        }
        const Gen_hepevt& thisGen = GenMgr(
            Panther_ID(mother.mother().get_ID()));
        p.relation().genHepevt(thisGen);
        dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(5);
        return;
      }
      return;
    }
  } else {
    std::cout << "[Zupanc] setMCtruthPi0: "
              << "Neither of two photons has link to gen_hepevt!"
              << std::endl;
  }
  return;
}


void setMCtruth(Particle &p)
{
  if (! &p.userInfo()) setUserInfo(p); // create UserInfo, if it does not exist
  // std::cout << " *** MC TRUTH *** " << std::endl;
  // std::cout << " Particle = " << p.lund() << " with " << p.nChildren()
  //           << " children." << std::endl;
  // std::cout << " Status before: Checked = "
  //           << dynamic_cast<UserInfo&>(p.userInfo()).genHepevtChecked()
  //           << " with Link = "
  //           << dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink()
  //           << std::endl;

  if (dynamic_cast<UserInfo&>(p.userInfo()).genHepevtChecked() == 1) {
    // std::cout << "Nothing to do - particle already checked!" << std::endl;
    return;
  }

  if (p.lund() == 111) {
    // std::cout << " Checking Pi0 " << std::endl;
    setMCtruthPi0(p);
    // std::cout << " Pi0 link ("
    //           << dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink()
    //           << ") " << std::endl;
    return;
  }

  // search for gen_hepevt link of FSP directly
  if (p.mdstCharged()) {
    //std::cout << " -> Charged track: " << std::endl;
    const Gen_hepevt & hep(get_hepevt(p.mdstCharged()));
    dynamic_cast<UserInfo&>(p.userInfo()).genHepevtChecked(1);
    if (hep) {
      p.relation().genHepevt(hep);
      if (p.lund() == hep.idhep()) {
        dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(1);
        // std::cout << " link found ("
        //           << dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink()
        //           << ") " << std::endl;
      } else {
        dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(2);
        // std::cout << " link found ("
        //           << dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink()
        //           << ") " << std::endl;
      }
      return;
    }
  } else if (p.mdstGamma()) {
    // std::cout << " -> Gamma: " << std::endl;
    const Gen_hepevt & hep(get_hepevt(p.mdstGamma()));
    dynamic_cast<UserInfo&>(p.userInfo()).genHepevtChecked(1);
    if (hep) {
      p.relation().genHepevt(hep);
      dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(1);
      // std::cout << " link found ("
      //           << dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink()
      //           << ") " << std::endl;
      return;
    }
  }

  // it is not mdstCharged or mdstGamma (combined particle)
  int nChildren = p.relation().nChildren();
  if (nChildren < 2)
    return;

  // special treatment for pi0
  // in case that one of the gammas doesn't have link to gen_hepevt but the
  // other one has set link for pi0 to the mother of that gamma.

  // set checked flag to 1
  dynamic_cast<UserInfo&>(p.userInfo()).genHepevtChecked(1);

  // check that all child particles have MC truth
  // std::cout << " Particle != FSP : nChidlren = "
  //           << p.relation().nChildren() << std::endl;
  for(int i=0; i<nChildren; ++i) {
    // std::cout << " [" << i << "]" << "lund = "
    //           << p.relation().child(i).lund() << " checked/link = "
    //           << dynamic_cast<UserInfo&>(p.child(i).userInfo()).genHepevtChecked()
    //           << "/" << dynamic_cast<UserInfo&>(p.child(i).userInfo()).genHepevtLink()
    //           << std::endl;
    if (!p.relation().child(i).genHepevt()) setMCtruth(p.relation().child(i));
    // std::cout << " [" << i << "]" << "lund = "
    //           << p.relation().child(i).lund() << " checked/link = "
    //           << dynamic_cast<UserInfo&>(p.child(i).userInfo()).genHepevtChecked()
    //           << "/" << dynamic_cast<UserInfo&>(p.child(i).userInfo()).genHepevtLink()
    //           << std::endl;
    if (!p.relation().child(i).genHepevt()) {
      // std::cout << " -> No link for child " << i << std::endl;
      return;
    }
    if (p.relation().child(i).genHepevt().idhep() == 0) {
      //std::cout << " -> Link = 0 for child " << i << std::endl;
      return;
    }
  }

  // Check that there are no clones
  for(int i=0; i<nChildren; ++i) {
    for(int j=i+1; j<nChildren; ++j) {
      // std::cout << "[" << i << "/" << j << "] = "
      //           << p.relation().child(i).genHepevt().get_ID()
      //           << " ?? " << p.relation().child(j).genHepevt().get_ID()
      //           << std::endl;
      if  (p.relation().child(i).genHepevt().get_ID() ==
          p.relation().child(j).genHepevt().get_ID()) {
        // std::cout << "Clones found! Exit." << std::endl;
        return;
      }
    }
  }

  // Find common mother, if exists (p has at least 2 children).
  std::vector<int> thisMothers;
  std::vector<int> otherMothers;

  fillMothers(p, thisMothers, otherMothers);
  if (thisMothers.size() == 0 || otherMothers.size() == 0) {
    return;
  }
  int motherID =  findCommonMother(p.nChildren(), thisMothers, otherMothers);

  if (!motherID) {
    std::cout << "[Zupanc] setMCtruth: "
              << "requesting particle with ID = 0! Exiting." << std::endl;
    return;
  }

  Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
  const Gen_hepevt& thisGen = GenMgr(Panther_ID(motherID));
  //std::cout << "Common mother: ID = " << motherID << " idhep = "
  //          << thisGen.idhep() << std::endl;

  // Set the relation.
  p.relation().genHepevt(thisGen);

  // Number and type of final state particles should be the same.
  // ID = 1; final state is correctly reconstructed (with particle ID)
  // ID = 2; one or more FSP are misidentified, but have common mother
  // ID = 3; FSP have common mother, but at least one
  //    massive particle is missing
  // ID = 4; FSP have common mother, but at least one
  //    massless particle is missing (FSR photon or neutrino)
  // ID = 5; final state includes pi0, which has only one daughter gamma with
  //    link to gen_hepevt
  // ID = 6; ID = 2 and 5 are true
  // ID = -1; common mother is virtual gamma (10022) or
  //    Upsilon(4S,5S) (300553,9000553)
  int motherIDhep = thisGen.idhep();
  if (motherIDhep == 10022 || motherIDhep == 300553 || motherIDhep == 9000553) {
    dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(-1);
    return;
  }

  // Ks doesn't have daughters in gen_hepevt table
  if (motherIDhep == 310) {
    dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(1);
    return;
  }

  std::vector<Particle> reconstructed;
  std::vector<Gen_hepevt> generated;

  // std::cout << "Appending Rec/Gen FSP" << std::endl;
  appendRecFSP(p, reconstructed);
  int a = appendGenFSP(thisGen, generated);
  if (a == 0) {
    dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(-2);
    return;
  }

  int truth = compareFinalStates(reconstructed, generated);
  dynamic_cast<UserInfo&>(p.userInfo()).genHepevtLink(truth);
}



void setMCtruth(std::vector<Particle> &plist) {
  if (plist.size() == 0) return;

  for(std::vector<Particle>::iterator i = plist.begin();
      i != plist.end(); ++i) {
    setMCtruth(*i);
  }
}

////////////////////////////////////////////////////////////////////////////////
// get decay chain for final stat particle

void genDecayChain(Particle p, int* dChain) {
  for(int i=0; i<=8; i++) dChain[i] = -1;

  if (p.relation().genHepevt()) {
    Gen_hepevt igen = p.relation().genHepevt();
    dChain[0] = igen.idhep();

    Gen_hepevt imot = igen.mother();
    if (imot) {
      dChain[1] = imot.idhep();
      dChain[2] = imot.daLast()-imot.daFirst()+1;

      Gen_hepevt immot = imot.mother();
      if (immot) {
        dChain[3] = immot.idhep();
        dChain[4] = immot.daLast()-immot.daFirst()+1;

        Gen_hepevt rg_mmmot = immot.mother();
        if (rg_mmmot) {
          dChain[5] = rg_mmmot.idhep();
          dChain[6] = rg_mmmot.daLast()-rg_mmmot.daFirst()+1;

          Gen_hepevt mrg_mmmot = rg_mmmot.mother();
          if (mrg_mmmot) {
            dChain[7] = mrg_mmmot.idhep();
            dChain[8] = mrg_mmmot.daLast()-mrg_mmmot.daFirst()+1;
          }
        }
      }
    }
  }
}

int IDhep(Particle &particle)
{
  if (! particle.genHepevt()) return 0;
  return particle.genHepevt().idhep();
}

int IDmom(Particle &particle)
{
  if (particle.relation().genHepevt()) {
    if (particle.relation().genHepevt().mother()) {
      return particle.relation().genHepevt().mother().idhep();
    }
  }
  return 0;
}

int NdecayProd(Particle &part)
{
  if (! part.genHepevt()) return 0;
  return part.genHepevt().daLast() - part.genHepevt().daFirst() +1;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
