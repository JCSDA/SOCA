/*
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"
#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"
#include "ufo/instantiateObsFilterFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory<soca::Traits>();
  oops::HofX3D<soca::Traits> hofx;
  run.execute(hofx);
  return 0;
}
