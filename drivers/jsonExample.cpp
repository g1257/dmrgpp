//-*-C++-*-

/** \ingroup main.JsonParser */
/*@{*/

/*! \file jsonParser.cpp 
 *
 *  
 *
 */

#include "JSON/JsonReader.h"
#include "JSON/JsonParser/DefaultContext.h"
#include <wchar.h>
#include <iostream>
#include <fstream>
#include "String.h"

int main(int argc,char *argv[]) {
  
  if (argc < 2) {
    std::cout << "Usage: "<<argv[0]<<" inputFileName\n";
    return -1;
  }

  PsimagLite::String    fileName(argv[1]);
/*  std::wifstream file;
  
  file.open(fileName.c_str());*/
  
//   JsonParser::JsonParser<JsonParser::DefaultContext> parser;
//   parser.filename = fileName;		
//   while(parser.parseChar(file));
	dca::JsonReader reader(fileName);
	int totalBins =0;
	totalBins <= reader["programSpecific"]["DCA"]["control"]["totalBins"];
	std::wcout <<" totalBins = "<<totalBins<<"\n";

	int numPointInBath = 0;
	numPointInBath <= reader.searchFor("numPointsInBath");
	std::wcout <<" numPointsInBath = "<<numPointInBath<<"\n";
  
  //std::wcout << parser.ctx;

  return 0;
}


/*@}*/

