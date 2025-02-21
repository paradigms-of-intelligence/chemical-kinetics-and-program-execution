/*
Copyright 2025 Google LLC

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <stdlib.h>
#include <gambit.h>

typedef ___setup_params_struct GAMBIT;

/* This must align with the symbol that got defined
   in the <name>_.c created by gsc --link. */
#define SCHEME_LIBRARY_LINKER ___LNK_tapes__py__interface__


___BEGIN_C_LINKAGE
extern ___mod_or_lnk SCHEME_LIBRARY_LINKER (___global_state);
___END_C_LINKAGE

GAMBIT* setup_gambit(){
  GAMBIT* g = malloc(sizeof(GAMBIT));
    
  ___setup_params_reset(g);
  g->version = ___VERSION;
  g->linker = SCHEME_LIBRARY_LINKER;

  ___setup(g);
  return g;
}

void cleanup_gambit(GAMBIT* g){
  ___cleanup();
  free(g);
}
