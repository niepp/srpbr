#ifndef __MODEL_H__
#define __MODEL_H__

#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

typedef std::vector<model_vertex_t> model_vertex_vec_t;
typedef std::vector<interp_vertex_t> interp_vertex_vec_t;
typedef std::vector<uint16_t> index_vec_t;

struct model_t 
{
	model_vertex_vec_t m_model_vertex;
	index_vec_t m_model_indices;
	interp_vertex_vec_t m_vertex_post;

	struct index_t
	{
		int iv, it, in; // index of vert / uv / normal
		index_t() : iv(-1), it(-1), in(-1) {}

		bool valid() const {
			return iv >= 0 && it >= 0 && in >= 0;
		}

		friend bool operator<(const index_t& l, const index_t& r) {
			assert(l.valid() && r.valid());
			return (l.iv < r.iv) ? true : (
				l.iv > r.iv ? false : (
					l.it < r.it ? true : (
						l.it > r.it ? false : l.in < r.in
						)
					)
				);
		}
	};

	struct triangle_t {
		index_t a, b, c;
	};

	model_t()
	{
	}

	int load(const std::string &mpath)
	{
		std::ifstream in;
		in.open(mpath, std::ifstream::in);
		if (in.fail()) {
			printf("load model failed\n");
			return -1;
		}

		auto slash_count = [](const std::string& instr, char flag) -> int {
			int scnt = 0;
			for (auto c : instr) {
				scnt += (c == flag);
			}
			return scnt;
		};

		auto read_point = [](std::istringstream& s, index_t &p) {
			char slash;
			s >> p.iv >> slash >> p.it >> slash >> p.in;
			// in wavefront obj all indices start at 1, not zero
			--p.iv;
			--p.it;
			--p.in;
			assert(p.valid());
		};

		std::vector<vector3_t> verts;
		std::vector<triangle_t> faces;  // index of vertex/uv/normal
		std::vector<vector3_t> norms;
		std::vector<texcoord_t> uvs;

		std::string line;
		while (!in.eof()) {
			std::getline(in, line);
			std::istringstream iss(line.c_str());
			char trash;
			if (!line.compare(0, 2, "v "))
			{
				iss >> trash;
				vector3_t v;
				iss >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}
			else if (!line.compare(0, 3, "vn "))
			{
				iss >> trash >> trash;
				vector3_t n;
				iss >> n.x >> n.y >> n.z;
				norms.push_back(n);
			}
			else if (!line.compare(0, 3, "vt "))
			{
				iss >> trash >> trash;
				texcoord_t tc;
				iss >> tc.u >> tc.v;
				uvs.push_back(tc);
			}
			else if (!line.compare(0, 2, "f "))
			{
				iss >> trash;
				int scnt = slash_count(line, '/');
				if (scnt != 6) {
					return -1;
				}
				triangle_t tri;
				read_point(iss, tri.a);
				read_point(iss, tri.b);
				read_point(iss, tri.c);
				faces.push_back(tri);
			}
		}

		std::cout << "# v# " << verts.size() << " f# " << faces.size() << " vt# " << uvs.size() << " vn# " << norms.size() << std::endl;

		int vidx = 0;
		std::map<index_t, int> vmap;
		for (auto& f : faces)
		{
			if (vmap.find(f.a) == vmap.end()) {
				vmap[f.a] = vidx++;
			}
			if (vmap.find(f.b) == vmap.end()) {
				vmap[f.b] = vidx++;
			}
			if (vmap.find(f.c) == vmap.end()) {
				vmap[f.c] = vidx++;
			}
		}

		m_model_vertex.resize(vmap.size());
		m_model_indices.resize(faces.size() * 3);
		m_vertex_post.resize(vmap.size());

		for (auto& kv : vmap)
		{
			model_vertex_t &vtx = m_model_vertex[kv.second];
			index_t idx = kv.first;
			vtx.pos = verts[idx.iv] * 0.01f;   // ue4 obj's unit is centimeter, convert to meter
			vtx.uv = uvs[idx.it];
			vtx.nor = norms[idx.in];
			vector3_t c = vtx.pos;
			c.normalize();
			vtx.color = c;
			std::swap(vtx.pos.y, vtx.pos.z);   // obj file format is right-hand coordinate system
			std::swap(vtx.nor.y, vtx.nor.z);   // obj file format is right-hand coordinate system
		}

		for (size_t i = 0; i < faces.size(); ++i)
		{
			int ia = vmap[faces[i].a];
			int ib = vmap[faces[i].b];
			int ic = vmap[faces[i].c];
			m_model_indices[i * 3 + 0] = ia;
			m_model_indices[i * 3 + 1] = ib;
			m_model_indices[i * 3 + 2] = ic;
		}

		return 0;

	}

};

#endif //__MODEL_H__


