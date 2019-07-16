/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "soca/Traits.h"
#include "soca/Transforms/instantiateBalanceOpFactory.h"
#include "oops/runs/EnsVariance.h"
#include "soca/Run/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  soca::instantiateBalanceOpFactory();
  oops::EnsVariance<soca::Traits> ensvariance;
  run.execute(ensvariance);
  return 0;
}