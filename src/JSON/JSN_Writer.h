//-*-C++-*-

/** \ingroup DCA */
/*@{*/

/*! \file JSN_Writer.h

*/

#ifndef dca_JSN_Writer_H
#define dca_JSN_Writer_H

namespace dca {

  class JSN {
  public:

    int width;
    bool doNothing;

    JSN(bool action=true):
      width (13),
      doNothing(action)
    {}
    
    std::string toString(const std::map<std::string, std::string>& map) {
      std::ostringstream buff;
      buff << "{";
      std::map<std::string,std::string>::const_iterator itr;
      for (itr = map.begin(); itr != map.end(); itr++) {
	if (itr != map.begin())
	  buff << ",";
	buff << " \"" << itr->first << "\": \"" << itr->second << "\"";
      }
      buff << "}";
      return buff.str();
    }

    template<typename T>
    std::string toString(const std::map<std::string, T>& map) {
      typedef std::map<std::string,T> MapType;
      typedef typename MapType::const_iterator CITR;

      std::ostringstream buff;
      buff << "{";
      for (CITR itr = map.begin(); itr != map.end(); itr++) {
	if (itr != map.begin())
	  buff << ",";
	buff << " \"" << itr->first << "\": " << itr->second;
      }
      buff << "}";
      return buff.str();
    }

    template<typename T>
    std::string toString(const std::map<std::string, std::vector<T> >& map) {

      typedef std::map<std::string,T> MapType;
      typedef typename MapType::const_iterator CITR;

      std::ostringstream buff;
      buff << "{";

      for (CITR itr = map.begin(); itr != map.end(); itr++) {
	if (itr != map.begin())
	  buff << ",";
	buff << " \"" << itr->first << "\": " <<  psimag::VectorLike::toString(itr->second);
      }
      buff << "}";
      return buff.str();
    }

    template<class T>
    inline
    void printVector(std::vector<T> vec,
		     std::string title,
		     std::ostream& os)  {
      os.precision(width);
      os << std::fixed; // scientific;
      os << "{ "
	 << " \'tile\': \'" << title   << "\', \n"
	 << " \'type\': \'Function1D', \n"
	 << " \'size\': " << vec.size() << ", \n";
      os << " \'data\': ";
      os << "array([[";
      for(size_t j=0; j<vec.size(); j++) 
	os << " " << std::setw(width) << vec[j] << ",\n";
      os << "]])\n";
      
      os << "}";
    }
 
    template<typename ObjType>
    void writeToFile(std::string fileName, 
		     std::string objName,
		     const ObjType& obj) {
      std::ofstream     file;
      std::string vtkFileName(fileName);
      vtkFileName += ".vts";
      file.open(fileName.c_str());
      obj.toJSN(file,objName,width);
      //obj.writeVTK(vtkFileName);
      file.close();
      
    }

    template<typename T>
    void writeToFile(std::string fileName, 
		     std::string objName,
		     const std::vector<T>& obj) {
      std::ofstream     file;
      file.open(fileName.c_str());
      printVector(obj,objName,file);
      file.close();

    }

    template<typename Type>
    void write(const std::vector<Type>& vec, 
	       std::string objName, 
	       std::string collectionName,
	       int         index,
	       std::string parentName) {
      if (doNothing) return;
      std::ostringstream fileNameBuff;
      fileNameBuff << "debug/" << collectionName << "/" << index << "/" << parentName << "/" << objName << ".py";
      writeToFile(fileNameBuff.str(),objName,vec);
    }
    


    template<typename ObjType>
    void write(const ObjType& obj, std::string fileName, std::string title="jsn object") {
      if (doNothing) return;
      writeToFile(fileName,title,obj);
    }
    
    template<typename ObjType>
    void write(const ObjType& obj, 
	       std::string objName, 
	       std::string collectionName,
	       int         index,
	       std::string parentName) {
      if (doNothing) return;
      std::ostringstream fileNameBuff;
      fileNameBuff << "debug/" << collectionName << "/" << index << "/" << parentName << "/" << objName << ".py";
      writeToFile(fileNameBuff.str(),objName,obj);
    }

    
    template<typename ObjType>
    void write(const ObjType& obj, 
	       std::string objName, 
	       std::string collectionName,
	       int         index,
	       std::string parentName,
	       std::string subCollectionName,
	       int         index2) {
      if (doNothing) return;
      std::ostringstream fileNameBuff;
      fileNameBuff << "debug/" << collectionName << "/" << index << "/" << parentName <<"/" 
		   << subCollectionName << "/" << index2 << "/" << objName << ".py";
      writeToFile(fileNameBuff.str(),objName,obj);
    }
    
    template<typename ObjType>
    void write(const ObjType& obj, 
	       std::string objName, 
	       std::string collectionName,
	       int         index,
	       std::string parentName,
	       std::string subCollectionName) {
      if (doNothing) return;
      std::ostringstream fileNameBuff;
      fileNameBuff << "debug/" << collectionName << "/" << index << "/" << parentName <<"/" 
		   << subCollectionName << "/" << objName<< ".py";
      std::string fileName(fileNameBuff.str());
      writeToFile(fileNameBuff.str(),objName,obj);
    }
    
  };
  
} // end namespace DCA


/*@}*/
#endif
