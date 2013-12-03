#include <iostream>
#include <vector>
#include <fstream>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

using namespace std;
int main(int argc, char **argv)
{
	MyMesh  mesh;
	// check command line options
	if (argc != 2) 
	{
		std::cerr << "Usage:  " << argv[0] << "save_normal infile\n";
		return 1;
	}
	// read mesh from stdin
	if ( ! OpenMesh::IO::read_mesh(mesh, argv[1]) )
	{
		std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
		return 1;
	}
	cout<<"input point size:"<<mesh.n_vertices()<<endl;
	cout<<"updating face normal..."<<endl;
	mesh.request_face_normals();
	mesh.update_face_normals();
	cout<<"updated face normal"<<endl;
	vector<OpenMesh::Vec3f> normals;
	int iter=0;
	for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
	{ 
		iter++;
		// do something with *v_it, v_it->, or v_it.handle()
		MyMesh::VertexFaceIter	vf_it;
		for (vf_it=mesh.vf_iter( v_it ); vf_it; ++vf_it)
		{
			//cout<<"iterating faces around the vertex"<<endl;
			OpenMesh::Vec3f n;
			n=mesh.normal(vf_it);
			if (n[0]==0 && n[1]==0 && n[2]==0)
			{
				cout<<"Warning: NULL normal!"<<endl;
				cout<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
			}
			n=n.normalize();
			normals.push_back(n);
			break;
		}
	}
	cout<<"iter:"<<iter<<endl;
	cout<<"added normal size:"<<normals.size()<<endl;
	cout<<"output"<<endl;
	ofstream outfile("standard_normal.xyzn");
	ofstream outfile1("tmp.xyzn");
	int index=0;
	for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
	{ 
		OpenMesh::Vec3f p=mesh.point(v_it);
		OpenMesh::Vec3f n=normals[index];
		outfile<<p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
		outfile1<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
		++index;
	}
	outfile.close();
	outfile1.close();
	return 0;
}