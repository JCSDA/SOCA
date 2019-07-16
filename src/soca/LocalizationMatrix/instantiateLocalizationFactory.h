/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SOCA_LOCALIZATIONMATRIX_INSTANTIATELOCALIZATIONFACTORY_H_
#define SOCA_LOCALIZATIONMATRIX_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/interface/LocalizationBase.h"
#include "soca/LocalizationMatrix/LocalizationMatrix.h"
#include "soca/Traits.h"

namespace soca {

void instantiateLocalizationFactory() {
}
// static oops::LocalizationMaker<soca::Traits,
//                                LocalizationMatrix> maker_("SOCA");
}  // namespace soca

#endif  // SOCA_LOCALIZATIONMATRIX_INSTANTIATELOCALIZATIONFACTORY_H_