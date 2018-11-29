#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
  for (Atom *atomi : m_atoms){
    if (atomi->position.x()<-m_systemSize.x()*0.5){
      atomi->position.components[0]+=m_systemSize.x();
      }
    if (atomi->position.x()>= m_systemSize.x()*0.5){
      atomi->position.components[0]-=m_systemSize.x();
	}
    
    if (atomi->position.y()<-m_systemSize.y()*0.5){
      atomi->position.components[1]+=m_systemSize.y();
    }
    if (atomi->position.y()>= m_systemSize.y()*0.5){
      atomi->position.components[1]-=m_systemSize.y();
    }
       
    if (atomi->position.z()<-m_systemSize.z()*0.5){
      atomi->position.components[2]+=m_systemSize.z();
    }
    if (atomi->position.z()>= m_systemSize.z()*0.5){
      atomi->position.components[2]-=m_systemSize.z();
    }
    
    for(Atom *atomj : m_atoms){
      if(atomi!=atomj){
      double dx=atomj->position.x()-atomi->position.x();
      if (dx>m_systemSize.x()*0.5){
	dx-=m_systemSize.x();
      }
      if (dx<= -m_systemSize.x()*0.5){
	dx+=m_systemSize.x();
      }
      double dy=atomj->position.y()-atomi->position.y();
      if (dy>m_systemSize.y()*0.5){
	dy-=m_systemSize.y();
      }
      if (dy<= -m_systemSize.y()*0.5){
	dy+=m_systemSize.y();
      }
      double dz=atomj->position.z()-atomi->position.z();
      if (dz>m_systemSize.z()*0.5){
	dz-=m_systemSize.z();
      }
      if (dz<= -m_systemSize.z()*0.5){
	dz+=m_systemSize.z();
      }
      } 
    }
  }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
  for (Atom*atomi:m_atoms){
    double px=atomi->mass()*atomi->velocity.x();
    double py=atomi->mass()*atomi->velocity.y();
    double pz=atomi->mass()*atomi->velocity.z();
    momentum.push_back((px, py, pz));
  }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

    for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
