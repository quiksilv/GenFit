/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef genfit_myProlateSpacepointMeasurement_h
#define genfit_myProlateSpacepointMeasurement_h

#include "ProlateSpacepointMeasurement.h"
#include "TrackCandHit.h"
#include "mySpacepointDetectorHit.h"


namespace genfit {

/** @brief Example class for a spacepoint measurement which can be created
 * from mySpacepointDetectorHit via the MeasurementFactory.
 *
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */
class myProlateSpacepointMeasurement : public ProlateSpacepointMeasurement {

 public:

  /** Default constructor for ROOT IO. */
  myProlateSpacepointMeasurement() :
     ProlateSpacepointMeasurement() {;}

  myProlateSpacepointMeasurement(const mySpacepointDetectorHit* detHit, const TrackCandHit* hit) :
    ProlateSpacepointMeasurement()
  {
    rawHitCoords_(0) = detHit->getPos()(0);
    rawHitCoords_(1) = detHit->getPos()(1);
    rawHitCoords_(2) = detHit->getPos()(2);
    rawHitCov_ = detHit->getCov();
    detId_ = hit->getDetId();
    hitId_ = hit->getHitId();


    largestErrorDirection_ = TVector3(1,0,0);
    ////largestErrorDirection_ = TVector3(0., 0. ,1.);

    this -> initG();
  }

  virtual myProlateSpacepointMeasurement* clone() const {return new myProlateSpacepointMeasurement(*this);}

  ClassDef(myProlateSpacepointMeasurement,1)
};
/** @} */

} /* End of namespace genfit */

#endif // genfit_myProlateSpacepointMeasurement_h
