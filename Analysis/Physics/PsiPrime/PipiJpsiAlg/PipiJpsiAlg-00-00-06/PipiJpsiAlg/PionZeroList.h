#ifndef Physics_Analysis_PionZeroList_H
#define Physics_Analysis_PionZeroList_H

#include <vector>
#include <algorithm>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
using namespace CLHEP;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
typedef std::vector<Hep3Vector> Vp3;
typedef std::vector<HepLorentzVector> Vp4;

class PionZero{
 public:
  PionZero(HepLorentzVector &vp_gam_1, int &index_gam_1, HepLorentzVector &vp_gam_2, int &index_gam_2);
  PionZero() {};
  int get_index(int i){ return i>0 ? m_index_2 : m_index_1; }
  double get_mass() { return m_inv_mass; }
  double get_goodness() { return m_goodness; }
  void cal_goodness(int good_method); // 0: mass diff; 1: chi^2; 2: ...???
  HepLorentzVector get_gam_vp(int i){return i>0 ? m_vp_2 : m_vp_1;}
  HepLorentzVector get_pi0_vp(){return m_vp_total;}
 private:
  int m_index_1, m_index_2, m_good_method;
  HepLorentzVector m_vp_1, m_vp_2, m_vp_total;
  double m_goodness; // desribe how good of a pi0
  double m_inv_mass;
};


class PionZeroList{
 public:
  PionZeroList(Vp4 &input_vp); // initial photon four-momentum list input
  void set_cut(int index, double cut);
  void sort(); // sort all pion0 list with goodness( lower is better )
  void reduce(); // delete all the pion0 with re-used photons
  void print(); // print all information in pion0 list
  void refresh(); // after reduce and sort, we should refresh the mom list
  int get_num_pi0(){ return m_num_pi0; }
  Vp4 get_pi0_list(){ return m_pi0_vp; }

  class f_less: public std::binary_function<PionZero, PionZero, bool>{
  public:
    result_type operator()(first_argument_type p1, second_argument_type p2){
      return (result_type) (p1.get_goodness() < p2.get_goodness());
    }
  };

 private:
  Vp4 m_gam_vp, m_pi0_vp;
  std::vector<PionZero> Vpi0_list;
  double m_low_cut, m_high_cut; // mass cuts of pion0
  int m_num_pi0;
};

#endif
