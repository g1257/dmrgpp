//-*-C++-*-

/** \ingroup main.JsonParser */
/*@{*/

/*! \file main.c 
 *
 *  
 *
 */

#ifndef  JsonParser_MAIN_H
#define  JsonParser_MAIN_H

#include "JsonParser.h"
#include "DefaultContext.h"
#include <wchar.h>
#include <iostream>
#include <fstream>
#include "String.h"

int main(int argc,char *argv[]) {
  
  if (argc < 2) {
    std::cout << "Usage: "<<argv[0]<<" inputFileName\n";
    return -1;
  }

  String    fileName(argv[1]);
  std::wifstream file;
  
  file.open(fileName.c_str());
  
  JsonParser::JsonParser<JsonParser::DefaultContext> parser;
  parser.filename = fileName;		
  while(parser.parseChar(file));
  
  std::wcout << parser.ctx;
  
  return 0;
}


/*@}*/
#endif
