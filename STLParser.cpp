#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <streambuf>

#include <string>
#include <vector>

#include <map>


float eps_3 = (float)0.001; // 10e-3
bool IsRelativelyEqual(const float& eps, const float& a, const float& b)
{
  float left = std::fabs(a - b);
  float right = (eps * (1 + (std::fabs(a) + fabs(b)) / 2.0));
  return (left <= right);
};

class Vector
{
  float m_x;
  float m_y;
  float m_z;
public:
  Vector() : m_x(0), m_y(0), m_z(0) {}
  Vector(float xp, float yp, float zp) : m_x(xp), m_y(yp), m_z(zp) {}
  float GetX() { return m_x; };
  float GetY() { return m_y; };
  float GetZ() { return m_z; };
  bool operator==(Vector other)
  {
    return (IsRelativelyEqual(eps_3, m_x, other.GetX()) &&
      IsRelativelyEqual(eps_3, m_y, other.GetY()) &&
      IsRelativelyEqual(eps_3, m_z, other.GetZ()));
  }

  // scalar product
  float operator * (Vector vec)
  {
    return (this->GetX() * vec.GetX()) + (this->GetY() * vec.GetY()) + (this->GetZ() * vec.GetZ());
  }

  //Norm
  Vector Norm()
  {
    float magnitude = sqrt((m_x * m_x) + (m_y * m_y) + (m_z * m_z));
    return Vector(m_x / magnitude, m_y / magnitude, m_z / magnitude);
  }

};

class Point {
  float m_x;
  float m_y;
  float m_z;

public:
  Point() : m_x(0), m_y(0), m_z(0) {}
  Point(float xp, float yp, float zp) : m_x(xp), m_y(yp), m_z(zp) {}
  float GetX() { return m_x; };
  float GetY() { return m_y; };
  float GetZ() { return m_z; };
  bool operator==(Point other)
  {
    return (m_x == other.GetX() && m_y == other.GetY() && m_z == other.GetZ());
  }

  Vector operator-(Point other)
  {
    return Vector(m_x - other.m_x, m_y - other.m_y, m_z - other.m_z);
  }
};

class Node {
  int m_index;
  float m_x;
  float m_y;
  float m_z;

public:
  Node() : m_index(-1), m_x(0), m_y(0), m_z(0) {}
  Node(int ind, float xp, float yp, float zp) : m_index(ind), m_x(xp), m_y(yp), m_z(zp) {}
  Point GetPoint() { return Point(m_x, m_y, m_z); };
  int GetIndex() { return m_index; }
};

class Triangle {
  int m_face_ind;
  int m_node_ind_1;
  int m_node_ind_2;
  int m_node_ind_3;
  Vector m_face_normal;

public:
  Triangle(int index, int n1, int n2, int n3, Vector normal) :
    m_face_ind(index), m_node_ind_1(n1), m_node_ind_2(n2), m_node_ind_3(n3), m_face_normal(normal) {}
  Vector GetNormal() { return m_face_normal; };
  int GetNode1() { return m_node_ind_1; };
  int GetNode2() { return m_node_ind_2; };
  int GetNode3() { return m_node_ind_3; };
};

struct stl_data {
  std::string name;
  std::vector<Triangle> triangles;

  stl_data(std::string namep) : name(namep) {}
};


class STL_Body
{
  std::vector<Node*> m_Geom_Nodes;
  std::vector<Triangle*>	m_Geom_Shapes;

public:
  STL_Body(std::vector<Node*> nodes, std::vector<Triangle*> shapes) : m_Geom_Nodes(nodes), m_Geom_Shapes(shapes) {};

  inline std::vector<Node*> GetNodes() { return m_Geom_Nodes; };
  inline std::vector<Triangle*> GetFaces() { return m_Geom_Shapes; };

};

std::vector<Point> m_vec_Vertices;
std::vector<Node*> m_Geom_Nodes;
std::map<int, Node*> m_Geom_Nodes_map;
std::vector<Triangle*>	m_Geom_Shapes;
std::vector<STL_Body> m_solid;

void ParseSTLFile(const std::string& stl_path);
void ParseAsciiFile(const std::string& stl_path);
void ParseBinaryFile(const std::string& stl_path);
float parse_float(std::ifstream& s);
Point parse_point(std::ifstream& s);
inline int GetNodesSize() { return m_Geom_Nodes.size(); };
bool isNumber(std::string str);
int AddAndGetPostionIndex(Point pt);
int AddTriangleToFacesList(int& faceIndex, int idx_1, int idx_2, int idx_3, Vector normal);
void GetTokensFromString(std::string inputString, std::vector<std::string>& tokens, std::string delimeter);

std::ostream& operator<<(std::ostream& out, Point p) {
  out << "v " << p.GetX() << " " << p.GetY() << " " << p.GetZ() << std::endl;
  return out;
}

//std::ostream& operator<<(std::ostream& out, Triangle t) {
//  out << "---- TRIANGLE ----" << std::endl;
//  out << t.GetNormal() << std::endl;
//  out << m_Geom_Nodes[t.GetNode1()] << std::endl;
//  out << m_Geom_Nodes[t.GetNode2()] << std::endl;
//  out << m_Geom_Nodes[t.GetNode3()] << std::endl;
//  return out;
//}

std::ostream& operator<<(std::ostream& out, Triangle t) {
  //out << "---- TRIANGLE ----" << std::endl;
  out << "f " << t.GetNode1() << ", " << t.GetNode2() << ", " << t.GetNode3() << std::endl;
  return out;
}

void write_output(std::string output)
{
  std::ofstream outdata; // outdata is like cin
  outdata.open(output);

  for (int i = 0; i < m_solid.size(); i++)
  {
    STL_Body body = m_solid[i];

    std::vector<Node*> nodes = body.GetNodes();
    std::vector<Triangle*> faces = body.GetFaces();
    outdata << "solid " << "  \n";

    for (int j = 0; j < nodes.size(); j++)
    {
      outdata << "v " << nodes[j]->GetPoint().GetX() << " " << nodes[j]->GetPoint().GetY() << " " << nodes[j]->GetPoint().GetZ() << " " << "\n";
    }
    for (int k = 0; k < faces.size(); k++)
    {
      outdata << "f " << faces[k]->GetNode1() << " " << faces[k]->GetNode2() << " " << faces[k]->GetNode3() << std::endl;
    }

    outdata << "endsolid \n";
  }
}

void write_Depth_for_facet(std::string output)
{
  std::ofstream outdata; // outdata is like cin
  outdata.open(output);

  for (int i = 0; i < m_solid.size(); i++)
  {
    STL_Body body = m_solid[i];

    std::vector<Node*> nodes = body.GetNodes();
    std::vector<Triangle*> faces = body.GetFaces();
    outdata << "solid " << "  \n";

    for (int k = 0; k < faces.size(); k++)
    {
      //outdata << "f " << faces[k]->GetNode1() << " " << faces[k]->GetNode2() << " " << faces[k]->GetNode3() << std::endl;
      int N1 = faces[k]->GetNode1();
      int N2 = faces[k]->GetNode2();
      int N3 = faces[k]->GetNode3();

      float depth = 0;
      Vector n = faces[k]->GetNormal().Norm();
      // Get any one point of the triangle
      Point TPt;
      //for (int j = 0; j < nodes.size(); j++)
      //{
      //  Node* Node = nodes[j];
      //  int id = Node->GetIndex();

      //  if (id == N1 || id == N2 || id == N3)
      //  {
      //    TPt = Node->GetPoint();
      //    break;
      //  }
      //}
      TPt = m_Geom_Nodes_map[N1]->GetPoint();

      //----
      for (int j = 0; j < nodes.size(); j++)
      {
        Node* Node = nodes[j];
        int id = Node->GetIndex();

        if (id == N1 || id == N2 || id == N3)
          continue;

        Point NPt = Node->GetPoint();
        Vector vTptNpt = NPt - TPt;
        float projection = vTptNpt * n;
        // depth in inward direction 
        if (projection <= 0.0)
        {
          float proj = std::fabs(projection);
          if (depth < proj)
            depth = proj;
        }

      }
      outdata << "f " << faces[k]->GetNode1() << " " << faces[k]->GetNode2() << " " << faces[k]->GetNode3() << " depth " << depth << std::endl;
    }

    outdata << "endsolid \n";
  }
}

int main()
{
  // Input File
  //std::string stl_file_name = "Bottle Holder 1 v2_bin.stl";
  std::string stl_file_name = "Fish tray T0009266 download free stl files 3d model for CNC wood carving.stl";
  //std::string stl_file_name = "Cube.stl";
  //-------------------------------------------------------------------------------------------------------------

  std::cout << "Reading STL data from file: " << stl_file_name << " ..." << std::endl;
  ParseSTLFile(stl_file_name);

  //-------------------------------------------------------------------------------------------------------------
  // Output File
  //std::string output_file_name = "triangleDepth_Bottle Holder 1 v2_bin.txt";
  std::string output_file_name = "triangleDepth_Fish tray T0009266 download free stl files 3d model for CNC wood carving.txt";
  //std::string output_file_name = "triangleDepth_Cube.txt";
  //-------------------------------------------------------------------------------------------------------------

  std::cout << "Writing STL data to file: " << output_file_name << " ..." << std::endl;
  //write_output(output_file_name);
  write_Depth_for_facet(output_file_name);
  return 0;
}

float parse_float(std::ifstream& s) {
  char f_buf[sizeof(float)];
  s.read(f_buf, 4);
  float* fptr = (float*)f_buf;
  return *fptr;
}
Vector parse_vector(std::ifstream& s) {
  float x = parse_float(s);
  float y = parse_float(s);
  float z = parse_float(s);
  return Vector(x, y, z);
}
Point parse_point(std::ifstream& s) {
  float x = parse_float(s);
  float y = parse_float(s);
  float z = parse_float(s);
  return Point(x, y, z);
}
void GetTokensFromString(std::string inputString, std::vector<std::string>& tokens, std::string delimeter)
{
  tokens.clear();
  std::size_t prev_pos = 0, pos;
  while ((pos = inputString.find_first_of(delimeter, prev_pos)) != std::string::npos)
  {
    if (pos > prev_pos)
      tokens.push_back(inputString.substr(prev_pos, pos - prev_pos));
    prev_pos = pos + 1;
  }
  if (prev_pos < inputString.length())
    tokens.push_back(inputString.substr(prev_pos, std::string::npos));
}
bool isNumber(std::string str)
{
  double d;
  std::istringstream is(str);
  is >> d;
  return !is.fail() && is.eof();
}
int AddAndGetPostionIndex(Point pt)
{
  int index = -1;
  auto iter = std::find(m_vec_Vertices.begin(), m_vec_Vertices.end(), pt);
  if (iter != m_vec_Vertices.end())
  {
    index = iter - m_vec_Vertices.begin();
  }
  else
  {
    m_vec_Vertices.push_back(pt);
    index = m_vec_Vertices.size() - 1;
  }
  return index;
}
int AddTriangleToFacesList(int& faceIndex, int idx_1, int idx_2, int idx_3, Vector normal)
{
  if (idx_1 == idx_2 || idx_2 == idx_3 || idx_3 == idx_1)
    return 1; // two vertices out of three are same. Resultant will be a line

  Triangle* triangle = new Triangle(faceIndex, idx_1, idx_2, idx_3, normal);
  m_Geom_Shapes.push_back(triangle);
  faceIndex++;
  return 0;
}
void ParseSTLFile(const std::string& stl_path)
{
  std::ifstream stlFileStream;

  stlFileStream.open(stl_path.c_str());
  if (!stlFileStream.is_open())
    return;

  bool isASCIIFormat = false;
  std::string line;
  double totalLines = 0;
  if (std::getline(stlFileStream, line))
  {
    std::vector<std::string> tempstr;
    GetTokensFromString(line, tempstr, " ");
    if (tempstr[0] == "solid") // Solid Body Loop
    {
      stlFileStream.seekg(0);
      for (; std::getline(stlFileStream, line);)
      {
        totalLines++;
        GetTokensFromString(line, tempstr, " ");
        if (!tempstr.empty() && tempstr[0] == "facet")
        {
          isASCIIFormat = true;
          break;
        }
        //tempstr.clear();
      }
      stlFileStream.clear();
    }
  }
  stlFileStream.close();

  if (isASCIIFormat)
    ParseAsciiFile(stl_path);
  else
    ParseBinaryFile(stl_path);
}
void ParseAsciiFile(const std::string& stl_path)
{
  std::ifstream stlFileStream;
  stlFileStream.open(stl_path.c_str());
  double currentLine = 0;
  for (std::string line; std::getline(stlFileStream, line);)
  {
    std::vector<std::string> tempstr;
    GetTokensFromString(line, tempstr, " ");
    if (tempstr[0] == "solid") // Solid Body Loop
    {
      int faceIndex = 0;
      for (; std::getline(stlFileStream, line);)
      {
        //tempstr.clear();
        GetTokensFromString(line, tempstr, " ");
        if (tempstr.size() == 0)
          continue;
        if (tempstr[0] == "facet")
        {
          bool isValid_number = isNumber(tempstr[2]);
          if (isValid_number) isValid_number = isNumber(tempstr[3]);
          else continue;
          if (isValid_number) isValid_number = isNumber(tempstr[4]);
          else continue;
          if (!isValid_number) continue;

          float x = stof(tempstr[2]);
          float y = stof(tempstr[3]);
          float z = stof(tempstr[4]);
          Vector faceNormal = Vector(x, y, z);

          int vertexCnt = 0;
          Point vertPt1;
          Point vertPt2;
          Point vertPt3;
          for (; std::getline(stlFileStream, line);)
          {
            GetTokensFromString(line, tempstr, " ");
            if (tempstr.size() == 0)
              continue;
            if (tempstr[0] == "vertex")
            {
              float x = stof(tempstr[1]);
              float y = stof(tempstr[2]);
              float z = stof(tempstr[3]);

              switch (vertexCnt)
              {
              case 0:
                vertPt1 = Point(x, y, z);
                break;
              case 1:
                vertPt2 = Point(x, y, z);
                break;
              case 2:
                vertPt3 = Point(x, y, z);
                break;
              }
              vertexCnt++;
            }
            else if (tempstr[0] == "endfacet")
              break;
          }
          if (!isValid_number)
            continue;

          int maxNodeIdx = GetNodesSize() - 1;
          int index_Pos1 = AddAndGetPostionIndex(vertPt1);
          int index_Pos2 = AddAndGetPostionIndex(vertPt2);
          int index_Pos3 = AddAndGetPostionIndex(vertPt3);

          if (index_Pos1 > maxNodeIdx)
          {
            Node* pNode = new Node(index_Pos1, vertPt1.GetX(), vertPt1.GetY(), vertPt1.GetZ());
            m_Geom_Nodes.push_back(pNode);
            m_Geom_Nodes_map[index_Pos1] = pNode;
          }
          if (index_Pos2 > maxNodeIdx)
          {
            Node* pNode = new Node(index_Pos2, vertPt2.GetX(), vertPt2.GetY(), vertPt2.GetZ());
            m_Geom_Nodes.push_back(pNode);
            m_Geom_Nodes_map[index_Pos2] = pNode;
          }
          if (index_Pos3 > maxNodeIdx)
          {
            Node* pNode = new Node(index_Pos3, vertPt3.GetX(), vertPt3.GetY(), vertPt3.GetZ());
            m_Geom_Nodes.push_back(pNode);
            m_Geom_Nodes_map[index_Pos3] = pNode;
          }

          AddTriangleToFacesList(faceIndex, index_Pos1, index_Pos2, index_Pos3, faceNormal);
        }
        else if (tempstr[0] == "endsolid")
          break;
      }

      if (m_Geom_Shapes.size())
      {
        m_solid.push_back(STL_Body(m_Geom_Nodes, m_Geom_Shapes));
      }

      m_vec_Vertices.clear();
      m_Geom_Nodes.clear();
      m_Geom_Shapes.clear();
    }
  }
  stlFileStream.close();
}
void ParseBinaryFile(const std::string& stl_path)
{
  std::ifstream stlFileStream(stl_path, std::ios::in | std::ios::binary);
  if (!stlFileStream)
    return;

  char header[81];
  unsigned int n_facets;

  char header_info[80] = "";
  char n_triangles[4];
  stlFileStream.read(header_info, 80);
  stlFileStream.read(n_triangles, 4);
  unsigned int* r = (unsigned int*)n_triangles;
  n_facets = *r;

  int faceIndex = 0;
  for (size_t i = 0; i < n_facets; ++i)
  {
    Vector faceNormal = parse_vector(stlFileStream);
    Point vertPt1 = parse_point(stlFileStream);
    Point vertPt2 = parse_point(stlFileStream);
    Point vertPt3 = parse_point(stlFileStream);
    char dummy[2];
    stlFileStream.read(dummy, 2);

    int maxNodeIdx = GetNodesSize() - 1;
    int index_Pos1 = AddAndGetPostionIndex(vertPt1);
    int index_Pos2 = AddAndGetPostionIndex(vertPt2);
    int index_Pos3 = AddAndGetPostionIndex(vertPt3);

    if (index_Pos1 > maxNodeIdx)
    {
      //m_Geom_Nodes.push_back(new Node(index_Pos1, vertPt1.GetX(), vertPt1.GetY(), vertPt1.GetZ()));
      Node* pNode = new Node(index_Pos1, vertPt1.GetX(), vertPt1.GetY(), vertPt1.GetZ());
      m_Geom_Nodes.push_back(pNode);
      m_Geom_Nodes_map[index_Pos1] = pNode;
    }
    if (index_Pos2 > maxNodeIdx)
    {
      //m_Geom_Nodes.push_back(new Node(index_Pos2, vertPt2.GetX(), vertPt2.GetY(), vertPt2.GetZ()));
      Node* pNode = new Node(index_Pos2, vertPt2.GetX(), vertPt2.GetY(), vertPt2.GetZ());
      m_Geom_Nodes.push_back(pNode);
      m_Geom_Nodes_map[index_Pos2] = pNode;
    }
    if (index_Pos3 > maxNodeIdx)
    {
      //m_Geom_Nodes.push_back(new Node(index_Pos3, vertPt3.GetX(), vertPt3.GetY(), vertPt3.GetZ()));
      Node* pNode = new Node(index_Pos3, vertPt3.GetX(), vertPt3.GetY(), vertPt3.GetZ());
      m_Geom_Nodes.push_back(pNode);
      m_Geom_Nodes_map[index_Pos3] = pNode;
    }
    AddTriangleToFacesList(faceIndex, index_Pos1, index_Pos2, index_Pos3, faceNormal);
  }

  if (m_Geom_Shapes.size())
  {
    m_solid.push_back(STL_Body(m_Geom_Nodes, m_Geom_Shapes));
  }

  stlFileStream.close();
}
