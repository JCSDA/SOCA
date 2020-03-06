/*
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/Model.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::Model<soca::Traits> tests;
  run.execute(tests);
  return 0;
}
