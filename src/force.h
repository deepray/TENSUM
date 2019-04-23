#ifndef __FORCE_H__
#define __FORCE_H__


/*!
  Structure to map boundary face tags to the type of force
*/
struct ForceData
{
   std::string name;
   std::vector<int> face_type;
};

/*!
  Structure to map boundary face to the value of the force evaluation
*/
struct Force
{
   std::vector<int> face;
   Vector           value;
};

#endif
