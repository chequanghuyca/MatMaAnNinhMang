//Chế Quang Huy - 1812340

#include <windows.h>
#include <math.h>
#include <iostream>
#include <glut.h>
#include <gl.h>

#define PI 3.141592653589793
#define COLORNUM 14

using namespace std;

class ColorRGB {
public:
	float rColor, gColor, bColor;
	void set(float red, float green, float blue) {
		rColor = red;
		gColor = green;
		bColor = blue;
	}
	void set(ColorRGB& c) {
		rColor = c.rColor;
		gColor = c.gColor;
		bColor = c.bColor;
	}
	ColorRGB() {
		rColor = gColor = bColor = 0;
	}
	ColorRGB(float red, float green, float blue) {
		rColor = red;
		gColor = green;
		bColor = blue;
	}
};

ColorRGB pink(199.0 / 255.0, 21.0 / 255.0, 133.0 / 255);
ColorRGB blue(0.0, 0.0, 1.0);
ColorRGB yellow(1.0, 1.0, 0.0);
ColorRGB mint_red(0.0, 0.3, 0.5);
ColorRGB green(0.0, 204.0 / 255.0, 0.0);

class Point3D {
	public:
		float x, y, z;
		void set(float dx, float dy, float dz) {
			x = dx; 
			y = dy; 
			z = dz;
		}
		void set(Point3D& p) {
			x = p.x; 
			y = p.y; 
			z = p.z;
		}
		Point3D() {
			x = y = z = 0; 
		}
		Point3D(float dx, float dy, float dz) {
			x = dx; 
			y = dy; 
			z = dz;
		}
};

class VectorXYZ {
	public:
		float x, y, z;
		void set(float dx, float dy, float dz) {
			x = dx; 
			y = dy; 
			z = dz;
		}
		void set(VectorXYZ& v) {
			x = v.x; 
			y = v.y; 
			z = v.z;
		}
		void flip() {
			x = -x; 
			y = -y; 
			z = -z;
		}
		void normalize() {
			float temp = sqrt(x * x + y * y + z * z);
			x = x / temp;
			y = y / temp;
			z = z / temp;
		}
		VectorXYZ() {
			x = y = z = 0; 
		}
		VectorXYZ(float dx, float dy, float dz) {
			x = dx; 
			y = dy; 
			z = dz;
		}
		VectorXYZ(VectorXYZ& v) {
			x = v.x; 
			y = v.y; 
			z = v.z;
		}
		VectorXYZ cross(VectorXYZ b) {
			VectorXYZ c(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
			return c;
		}
		float dot(VectorXYZ b) {
			return x * b.x + y * b.y + z * b.z;
		}
};

class VertexID {
	public:
		int vertIndex;
		int colorIndex;
};

class Face {
	public:
		int 		nVerts;
		VertexID*	vert;
		VectorXYZ 	v_newell;
		int 		faceID;

		Face() {
			faceID = 0;
			nVerts = 0;
			vert = NULL;
		}
		~Face() {
			if (vert != NULL) {
				delete[] vert;
				vert = NULL;
			}
			nVerts = 0;
		}
};

class Mesh {
	public:
		int numVerts;
		Point3D* pt;
		int numFaces;
		Face* face;
		GLfloat diffWhite[4] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat diffMain[4] = { 0.0, 0.0, 0.0, 1.0 };

	public:
		Mesh() {
			numVerts = 0;
			pt = NULL;
			numFaces = 0;
			face = NULL;
		}
		~Mesh() {
			if (pt != NULL) {
				delete[] pt;
			}
			if (face != NULL) {
				delete[] face;
			}
			numVerts = 0;
			numFaces = 0;
		}

		void KhoiBanNguyet(int N, float H, float H1, float R, float R1, float R2, float fGocQuay, float fGocQuay1, float shift);
		void KhoiHopXoay(int N, float H, float H1, float R, float R1, float fGocQuay);
		void HinhTru(int N, float H, float R);
		void HinhHop(float fLength, float fHeight, float fWidth);
		void SetMau(int colorIndex);
		void KhoiTaoMau(ColorRGB& c);
		void DrawKhung();
		void ToMau();
		void Draw(bool mode, float amb[], float spec[], float shi);
		void TinhGocQuay();
};

double ArrayColor[COLORNUM][3] = { { 0.6117664,0.15294 , 0.690196 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 },
	{ 151.0 / 255.0,255.0 / 255 , 255.0 / 255 }, { 199.0 / 255.0,21.0 / 255 , 133.0 / 255 },
	{ 0 / 255.0,191.0 / 255 , 255.0 / 255 }, { 0 / 255.0,250.0 / 255 , 154.0 / 255 },
	{ 0.3, 0.3, 0.3 }, { 0 / 255.0,250.0 / 255 , 154.0 / 255 }, { 0.9,  0.9, 0.9 },
	{ 1.0, 0.5,  0.5 }, { 0.6117664,0.15294 , 0.690196 }
};


void Mesh::HinhTru(int N, float H, float R) {
	numVerts = 2 * N + 2;
	pt = new Point3D[numVerts];
	pt[2 * N].set(0, H / 2, 0);
	pt[2 * N + 1].set(0, -H / 2, 0);

	float gocQuay = 2 * PI / N;
	for (int i = 0; i < N; i++) {
		float x = R * sin(gocQuay * i);
		float z = -R * cos(gocQuay * i);
		pt[i].set(x, H / 2, z);
		pt[i + N].set(x, -H / 2, z);
	}

	numFaces = 3 * N;
	face = new Face[numFaces];
	for (int i = 0; i < numFaces; i++) {
		if (i < N - 1) {
			face[i].nVerts = 3;
			face[i].vert = new VertexID[face[i].nVerts];
			face[i].vert[0].vertIndex = 2 * N;
			face[i].vert[1].vertIndex = i + 1;
			face[i].vert[2].vertIndex = i;
		}
		else if (i == N - 1) {
			face[i].nVerts = 3;
			face[i].vert = new VertexID[face[i].nVerts];
			face[i].vert[0].vertIndex = 2 * N;
			face[i].vert[1].vertIndex = 0;
			face[i].vert[2].vertIndex = i;
		}
		else if (i < 2 * N - 1) {
			face[i].nVerts = 3;
			face[i].vert = new VertexID[face[i].nVerts];
			face[i].vert[0].vertIndex = i + 1;
			face[i].vert[1].vertIndex = 2 * N + 1;
			face[i].vert[2].vertIndex = i;
		}
		else if (i == 2 * N - 1) {
			face[i].nVerts = 3;
			face[i].vert = new VertexID[face[i].nVerts];
			face[i].vert[0].vertIndex = N;
			face[i].vert[1].vertIndex = 2 * N + 1;
			face[i].vert[2].vertIndex = i;
		}
		else if (i < 3 * N - 1) {
			face[i].nVerts = 4;
			face[i].vert = new VertexID[face[i].nVerts];
			face[i].vert[0].vertIndex = i - 2 * N + 1;
			face[i].vert[1].vertIndex = i - N + 1;
			face[i].vert[2].vertIndex = i - N;
			face[i].vert[3].vertIndex = i - 2 * N;
		}
		else {
			face[i].nVerts = 4;
			face[i].vert = new VertexID[face[i].nVerts];
			face[i].vert[0].vertIndex = 0;
			face[i].vert[1].vertIndex = N;
			face[i].vert[2].vertIndex = 2 * N - 1;
			face[i].vert[3].vertIndex = N - 1;
		}
	}
}

void Mesh::KhoiBanNguyet(int N, float H, float H1, float R, float R1, float R2, float fGocQuay, float fGocQuay1, float shift) {
	int i;
	int index = 0;
	float x, y, z;
	float phi = (fGocQuay * PI / 180) / N;
	float teta = (fGocQuay1 * PI / 180) / N;
	numVerts = 6 * N + 12;
	pt = new Point3D[numVerts];

	float gocQuay = (-fGocQuay * 0.5) * PI / 180;
	for (i = 0; i <= N; i++) {
		x = R * cos(phi * i + gocQuay);
		z = R * sin(phi * i + gocQuay);
		y = 0;
		pt[i].set(x, y, z);
		y = H;
		pt[i + N + 1].set(x, y, z);
		x = R2 * cos(phi * i + gocQuay);
		z = R2 * sin(phi * i + gocQuay);
		y = 0;
		pt[i + 2 * (N + 1)].set(x, y, z);
		y = H1;
		pt[i + 3 * (N + 1)].set(x, y, z);
	}

	gocQuay = -(fGocQuay1 * 0.5) * PI / 180;
	for (i = 0; i <= N; i++) {
		x = R1 * cos(gocQuay + teta * i) + shift;
		z = R1 * sin(gocQuay + teta * i);
		y = H;
		pt[i + 4 * (N + 1)].set(x, y, z);
		y = H1;
		pt[i + 5 * (N + 1)].set(x, y, z);
	}

	x = 0.9 * R2 * cos(-fGocQuay * 0.5 * PI / 180);
	z = 0.9 * R2 * sin(-fGocQuay * 0.5 * PI / 180);
	pt[N + 1 + 5 * (N + 1)].set(x, H, z);
	pt[N + 2 + 5 * (N + 1)].set(x, H1, z);
	pt[N + 5 + 5 * (N + 1)].set(x, 0, z);
	x = 0.9 * R2 * cos(fGocQuay * 0.5 * PI / 180);
	z = 0.9 * R2 * sin(fGocQuay * 0.5 * PI / 180);
	pt[N + 3 + 5 * (N + 1)].set(x, H, z);
	pt[N + 4 + 5 * (N + 1)].set(x, H1, z);
	pt[N + 6 + 5 * (N + 1)].set(x, 0, z);
	numFaces = 6 * N + 10;
	face = new Face[numFaces];

	for (int i = 0; i < N; i++) {
		face[index].nVerts = 4;
		face[index].vert = new VertexID[face[index].nVerts];
		face[index].vert[0].vertIndex = i + N + 1;
		face[index].vert[1].vertIndex = i;
		face[index].vert[2].vertIndex = i + 1;
		face[index].vert[3].vertIndex = i + 2 + N;
		index++;
	}
	for (i = 0; i < N; i++) {
		face[index].nVerts = 4;
		face[index].vert = new VertexID[face[index].nVerts];
		face[index].vert[0].vertIndex = i + 3 * (N + 1);
		face[index].vert[1].vertIndex = i + 2 * (N + 1);
		face[index].vert[2].vertIndex = i + 2 * (N + 1) + 1;
		face[index].vert[3].vertIndex = i + 3 * (N + 1) + 1;
		index++;
	}
	for (i = 0; i < N; i++) {
		face[index].nVerts = 4;
		face[index].vert = new VertexID[face[index].nVerts];
		face[index].vert[0].vertIndex = i + 1;
		face[index].vert[1].vertIndex = i;
		face[index].vert[2].vertIndex = i + 2 * (N + 1);
		face[index].vert[3].vertIndex = i + 1 + 2 * (N + 1);
		index++;
	}
	for (i = 0; i < N; i++) {
		face[index].nVerts = 4;
		face[index].vert = new VertexID[face[index].nVerts];
		face[index].vert[0].vertIndex = i + N + 1 + 3 * (N + 1);
		face[index].vert[1].vertIndex = i + N + 1;
		face[index].vert[2].vertIndex = i + N + 2;
		face[index].vert[3].vertIndex = i + N + 2 + 3 * (N + 1);
		index++;
	}
	for (i = 0; i < N; i++) {
		face[index].nVerts = 4;
		face[index].vert = new VertexID[face[index].nVerts];
		face[index].faceID = 1;
		face[index].vert[0].vertIndex = i + 5 * (N + 1) + 1;
		face[index].vert[1].vertIndex = i + 5 * (N + 1);
		face[index].vert[2].vertIndex = i + 4 * (N + 1);
		face[index].vert[3].vertIndex = i + 4 * (N + 1) + 1;
		index++;
	}
	for (i = 0; i < N; i++) {
		face[index].nVerts = 4;
		face[index].vert = new VertexID[face[index].nVerts];
		face[index].vert[0].vertIndex = i + 3 * (N + 1) + 1;
		face[index].vert[1].vertIndex = i + 3 * (N + 1);
		face[index].vert[2].vertIndex = i + 5 * (N + 1);
		face[index].vert[3].vertIndex = i + 5 * (N + 1) + 1;
		index++;
	}

	face[index].nVerts = 4;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].vert[0].vertIndex = 4 * (N + 1);
	face[index].vert[1].vertIndex = 5 * (N + 1);
	face[index].vert[2].vertIndex = N + 2 + 5 * (N + 1);
	face[index].vert[3].vertIndex = N + 1 + 5 * (N + 1);
	index++;

	face[index].nVerts = 4;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].vert[0].vertIndex = 4 * (N + 1) + N;
	face[index].vert[1].vertIndex = N + 3 + 5 * (N + 1);
	face[index].vert[2].vertIndex = N + 4 + 5 * (N + 1);
	face[index].vert[3].vertIndex = 4 * (N + 1) + 2 * N + 1;
	index++;

	face[index].nVerts = 4;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].faceID = 1;
	face[index].vert[0].vertIndex = 0;
	face[index].vert[1].vertIndex = N + 1;
	face[index].vert[2].vertIndex = N + 1 + 5 * (N + 1);
	face[index].vert[3].vertIndex = N + 5 + 5 * (N + 1);
	index++;

	face[index].nVerts = 4;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].faceID = 1;
	face[index].vert[0].vertIndex = N;
	face[index].vert[1].vertIndex = N + N + 1;
	face[index].vert[2].vertIndex = N + 3 + 5 * (N + 1);
	face[index].vert[3].vertIndex = N + 6 + 5 * (N + 1);
	index++;

	face[index].nVerts = 4;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].faceID = 1;
	face[index].vert[0].vertIndex = N + 5 + 5 * (N + 1);
	face[index].vert[1].vertIndex = N + 2 + 5 * (N + 1);
	face[index].vert[2].vertIndex = 3 * (N + 1);
	face[index].vert[3].vertIndex = 2 * (N + 1);
	index++;

	face[index].nVerts = 4;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].faceID = 1;
	face[index].vert[0].vertIndex = 2 * (N + 1) + N;
	face[index].vert[1].vertIndex = N + 5 * (N + 1) + 6;
	face[index].vert[2].vertIndex = 3 + 6 * (N + 1);
	face[index].vert[3].vertIndex = 3 * (N + 1) + N;
	index++;

	face[index].nVerts = 3;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].vert[0].vertIndex = N + 1;
	face[index].vert[1].vertIndex = 4 * (N + 1);
	face[index].vert[2].vertIndex = N + 1 + 5 * (N + 1);
	index++;

	face[index].nVerts = 3;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].vert[0].vertIndex = 4 * (N + 1) + N;
	face[index].vert[1].vertIndex = N + 1 + N;
	face[index].vert[2].vertIndex = N + 3 + 5 * (N + 1);
	index++;

	face[index].nVerts = 3;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].vert[0].vertIndex = N + 5 * (N + 1) + 2;
	face[index].vert[1].vertIndex = 5 * (N + 1);
	face[index].vert[2].vertIndex = 3 * (N + 1);
	index++;

	face[index].nVerts = 3;
	face[index].vert = new VertexID[face[index].nVerts];
	face[index].vert[0].vertIndex = 5 * (N + 1) + N;
	face[index].vert[1].vertIndex = N + 4 + 5 * (N + 1);
	face[index].vert[2].vertIndex = 3 * (N + 1) + N;
	index++;
};

void Mesh::KhoiHopXoay(int N, float H, float H1, float R, float R1, float fGocQuay) {
	const float gocStep = fGocQuay / N;
	const float phiStep = gocStep * PI / 180;
	const float firstGocQuay = (-fGocQuay * 0.5) * PI / 180;
	numVerts = N * 4 + 6;
	pt = new Point3D[numVerts];
	pt[0].set(0, 0, 0);
	pt[1].set(0, H, 0);

	for (int i = 0; i <= N; i++) {
		const float gocQuay = firstGocQuay + i * phiStep;
		pt[i + 2].set(R1 * cos(gocQuay), 0, R1 * sin(gocQuay));
		pt[i + 2 + N + 1].set(R1 * cos(gocQuay), H1, R1 * sin(gocQuay));
		pt[i + 2 + 2 * (N + 1)].set(R * cos(gocQuay), H, R * sin(gocQuay));
		pt[i + 2 + 3 * (N + 1)].set(R * cos(gocQuay), H1, R * sin(gocQuay));
	}

	numFaces = N * 5 + 4;
	face = new Face[numFaces];
	int curIndex = 0;

	for (int i = 0; i < N; i++) {
		face[curIndex].nVerts = 4;
		face[curIndex].vert = new VertexID[face[curIndex].nVerts];
		face[curIndex].faceID = 1;
		face[curIndex].vert[0].vertIndex = i + 2;
		face[curIndex].vert[1].vertIndex = i + N + 3;
		face[curIndex].vert[2].vertIndex = i + N + 4;
		face[curIndex].vert[3].vertIndex = i + 3;
		curIndex++;

		face[curIndex].nVerts = 4;
		face[curIndex].vert = new VertexID[face[curIndex].nVerts];
		face[curIndex].vert[0].vertIndex = i + N + 4;
		face[curIndex].vert[1].vertIndex = i + N + 3;
		face[curIndex].vert[2].vertIndex = i + 2 + 3 * (N + 1);
		face[curIndex].vert[3].vertIndex = i + N + 4 + 2 * (N + 1);
		curIndex++;

		face[curIndex].nVerts = 3;
		face[curIndex].vert = new VertexID[face[curIndex].nVerts];
		face[curIndex].vert[0].vertIndex = i + 2 + 2 * (N + 1);
		face[curIndex].vert[1].vertIndex = 1;
		face[curIndex].vert[2].vertIndex = i + 3 + 2 * (N + 1);
		curIndex++;

		face[curIndex].nVerts = 4;
		face[curIndex].vert = new VertexID[face[curIndex].nVerts];
		face[curIndex].vert[0].vertIndex = i + 3 + 2 * (N + 1) + N;
		face[curIndex].vert[1].vertIndex = i + 2 + 2 * (N + 1);
		face[curIndex].vert[2].vertIndex = i + 3 + 2 * (N + 1);
		face[curIndex].vert[3].vertIndex = i + 4 + 2 * (N + 1) + N;
		curIndex++;
	}

	face[curIndex].nVerts = 4;
	face[curIndex].vert = new VertexID[face[curIndex].nVerts];
	face[curIndex].vert[0].vertIndex = 2 + N + 2 * (N + 1);
	face[curIndex].vert[1].vertIndex = 1;
	face[curIndex].vert[2].vertIndex = 0;
	face[curIndex].vert[3].vertIndex = 2 + N;
	curIndex++;

	face[curIndex].nVerts = 4;
	face[curIndex].vert = new VertexID[face[curIndex].nVerts];
	face[curIndex].vert[0].vertIndex = 2;
	face[curIndex].vert[1].vertIndex = 0;
	face[curIndex].vert[2].vertIndex = 1;
	face[curIndex].vert[3].vertIndex = 2 * (N + 1) + 2;
	curIndex++;

	face[curIndex].nVerts = 4;
	face[curIndex].vert = new VertexID[face[curIndex].nVerts];
	face[curIndex].vert[0].vertIndex = 2 + 2 * (N + 1) + N + 1;
	face[curIndex].vert[1].vertIndex = 3 + N;
	face[curIndex].vert[2].vertIndex = 2;
	face[curIndex].vert[3].vertIndex = 2 * (N + 1) + 2;
	curIndex++;

	face[curIndex].nVerts = 4;
	face[curIndex].vert = new VertexID[face[curIndex].nVerts];
	face[curIndex].vert[0].vertIndex = 2 * N + 2 + 1 + 2 * (N + 1);
	face[curIndex].vert[1].vertIndex = 2 + N + 2 * (N + 1);
	face[curIndex].vert[2].vertIndex = 2 + N;
	face[curIndex].vert[3].vertIndex = N + N + 3;
	curIndex++;
};

void Mesh::HinhHop(float fLength, float fHeight, float fWidth) {
	numVerts = 8;
	numFaces = 6;
	pt = new Point3D[numVerts];
	face = new Face[numFaces];

	Point3D vertices[8] = {
		Point3D(-fLength / 2, fHeight / 2, fWidth / 2),
		Point3D(fLength / 2, fHeight / 2, fWidth / 2),
		Point3D(fLength / 2, fHeight / 2, -fWidth / 2),
		Point3D(-fLength / 2, fHeight / 2, -fWidth / 2),
		Point3D(-fLength / 2, -fHeight / 2, fWidth / 2),
		Point3D(fLength / 2, -fHeight / 2, fWidth / 2),
		Point3D(fLength / 2, -fHeight / 2, -fWidth / 2),
		Point3D(-fLength / 2, -fHeight / 2, -fWidth / 2)
	};

	for (int i = 0; i < numVerts; i++) {
		pt[i] = vertices[i];
	}

	VertexID faces[6][4] = {
		{ { 3 },{ 7 },{ 4 },{ 0 } },
		{ { 3 },{ 2 },{ 6 },{ 7 } },
		{ { 1 },{ 5 },{ 6 },{ 2 } },
		{ { 0 },{ 4 },{ 5 },{ 1 } },
		{ { 6 },{ 5 },{ 4 },{ 7 } },
		{ { 2 },{ 3 },{ 0 },{ 1 } }
	};

	for (int i = 0; i < numFaces; i++) {
		face[i].nVerts = 4;
		face[i].vert = new VertexID[face[i].nVerts];
		for (int j = 0; j < face[i].nVerts; j++) {
			face[i].vert[j] = faces[i][j];
		}
	}
}

void Mesh::KhoiTaoMau(ColorRGB& c) {
	float diffuse[3] = { c.rColor, c.gColor, c.bColor };
	copy(diffuse, diffuse + 3, diffMain);
}

void Mesh::SetMau(int colorIdx) {
	for (int i = 0; i < numFaces; i++) {
		for (int v = 0; v < face[i].nVerts; v++) {
			face[i].vert[v].colorIndex = colorIdx;
		}
	}
}

void Mesh::DrawKhung() {
	glColor3f(0, 0, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int f = 0; f < numFaces; f++) {
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++) {
			int iv = face[f].vert[v].vertIndex;
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

void Mesh::TinhGocQuay() {
	float mx;
	float my;
	float mz;
	for (int f = 0; f < numFaces; f++) {
		mx = 0;
		my = 0;
		mz = 0;
		for (int i = 0; i < face[f].nVerts; i++) {
			int Index = face[f].vert[i].vertIndex;
			int IndexNext = face[f].vert[(i + 1) % face[f].nVerts].vertIndex;
			mx += (pt[Index].y - pt[IndexNext].y) * (pt[Index].z + pt[IndexNext].z);
			my += (pt[Index].z - pt[IndexNext].z) * (pt[Index].x + pt[IndexNext].x);
			mz += (pt[Index].x - pt[IndexNext].x) * (pt[Index].y + pt[IndexNext].y);
		}
		face[f].v_newell.set(mx, my, mz);
		face[f].v_newell.normalize();
	}
}

void Mesh::ToMau() {
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	for (int f = 0; f < numFaces; f++) {
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++) {
			int iv = face[f].vert[v].vertIndex;
			int ic = face[f].vert[v].colorIndex;
			glColor3f(ArrayColor[ic][0], ArrayColor[ic][1], ArrayColor[ic][2]);
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

void Mesh::Draw(bool mode, float amb[], float spec[], float shi) {
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	for (int f = 0; f < numFaces; f++) {
		(face[f].faceID == 1 && mode == true) ? glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffWhite) : glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffMain);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, amb);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shi);
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++) {
			int iv = face[f].vert[v].vertIndex;
			glNormal3f(face[f].v_newell.x, face[f].v_newell.y, face[f].v_newell.z);
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

float* MaTran(float* m1, float* m2, float res[16]) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			float sum = 0;
			for (int k = 0; k < 4; k++) {
				sum += m1[4 * i + k] * m2[4 * k + j];
			}
			res[4 * i + j] = sum;
		}
	}
	return res;
}

float camGocQuay = 60, camChieuRong = 30, camDis = 75;
bool view2D = false, isKhung = false, mode = false;
float xoay = 0, khoiBanNguyet_xoay = 0, PRE = 0, ANG = 0;
int	scrWidth = 1200, scrHeight = 600, k = 0;

VectorXYZ v_0(cos(-PI / 6), 0, sin(-PI / 6));
VectorXYZ v_xoay;
Mesh khoiBanNguyet, khoiHopXoay, hinhHop, hinhTru, base;
GLfloat amb[4] = { 0.3, 0.3, 0.3, 1.0 }, spec[4] = { 1.0, 1.0, 1.0, 1.0 }, shi = 90.5;

float modelview[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
float R[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
float v[16];

void control(unsigned char key, int x, int y) {
	switch (key) {
		case '1':
			xoay += 3;
			if (xoay > 360) {
				xoay -= 360;
			}
			break;
		case '2':
			xoay -= 3;
			if (xoay < -360) {
				xoay += 360;
			}
			break;
		case '3':
			mode = !mode;
			break;
		case '+':
			camDis += 0.8;
			break;
		case '-':
			camDis -= 0.8;
			break;
		case 'w':
		case 'W':
			isKhung = !isKhung;
			break;
		case 'v':
		case 'V':
			view2D = !view2D;
			break;
	}
	glutPostRedisplay();
}

void controlKeyboard(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_UP:
			camChieuRong += 1.0;
			break;
		case GLUT_KEY_DOWN:
			camChieuRong -= 1.0;
			break;
		case GLUT_KEY_RIGHT:
			camGocQuay -= 2.0;
			break;
		case GLUT_KEY_LEFT:
			camGocQuay += 2.0;
			break;
	}
	glutPostRedisplay();
}

void updateGocQuay() {
	if (xoay >= 360) {
		xoay -= 360;
	}
	if (xoay <= -360) {
		xoay += 360;
	}

	v_xoay.set(-3.6 * cos(xoay * PI / 180) + 7.2, 0, -3.6 * sin(xoay * PI / 180));
	v_xoay.normalize();

	float arr[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
	float gocQuay;

	switch ((int)xoay / 60) {
	case -6:
	case -5:
	case -4:
	case -3:
	case -2:
	case -1:
		if (xoay <= -60 && xoay >= -300) {
			ANG = acos(v_xoay.dot(v_0)) * 180 / PI;
			gocQuay = ANG - PRE;
			R[0] = cos(gocQuay * PI / 180);
			R[2] = -sin(gocQuay * PI / 180);
			R[8] = sin(gocQuay * PI / 180);
			R[10] = cos(gocQuay * PI / 180);
			MaTran(R, modelview, arr);
			for (int i = 0; i < 16; i++) {
				modelview[i] = arr[i];
			}
			PRE = ANG;
		}
		else {
			PRE = 0;
		}
		break;
	case 0:
		PRE = 0;
		break;
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
		if (xoay >= 60 && xoay <= 300) {
			ANG = acos(v_xoay.dot(v_0)) * 180 / PI;
			gocQuay = ANG - PRE;
			R[0] = cos(gocQuay * PI / 180);
			R[2] = -sin(gocQuay * PI / 180);
			R[8] = sin(gocQuay * PI / 180);
			R[10] = cos(gocQuay * PI / 180);
			MaTran(R, modelview, arr);
			for (int i = 0; i < 16; i++) {
				modelview[i] = arr[i];
			}
			PRE = ANG;
		}
		else {
			PRE = 60;
		}
		break;
	default:
		break;
	}
}

void drawNen(float alpha) {
	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColor4f(0.3, 1.0, 1.0, alpha);
	glBegin(GL_QUADS);
	glNormal3f(0.0f, 1.0f, 0.0f);
	for (float x = -31.4; x < 30; x += 2.1) {
		for (float z = -30; z < 30; z += 2.1) {
			glVertex3f(x, -2, z);
			glVertex3f(x + 2, -2, z);
			glVertex3f(x + 2, -2, z + 2);
			glVertex3f(x, -2, z + 2);
		}
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

void drawAllMesh() {
	float m[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
	updateGocQuay();

	if (!isKhung) {
		glPushMatrix();
		glTranslatef(7.2, 0, 0);
		glRotatef(xoay, 0, 1, 0);
		glTranslatef(-7.2, 0, 0);
		glTranslatef(3.6, 6.6, 0);
		hinhTru.KhoiTaoMau(green);
		hinhTru.Draw(mode, amb, spec, shi);
		glPopMatrix();
		glPushMatrix();
		glTranslatef(7.2, 5.8, 0);
		glRotatef(xoay, 0, 1, 0);
		khoiHopXoay.KhoiTaoMau(yellow);
		khoiHopXoay.Draw(mode, amb, spec, shi);
		glPopMatrix();
		glPushMatrix();
		glTranslatef(7.2, 2.2 + 1.5 * 0.5 + 3.6, 0);
		glRotatef(xoay, 0, 1, 0);
		hinhHop.KhoiTaoMau(yellow);
		hinhHop.Draw(mode, amb, spec, shi);
		glPopMatrix();
		glPushMatrix();
		glMultMatrixf(modelview);
		glPushMatrix();
		glTranslatef(0, 3.6 * 0.5, 0);
		base.KhoiTaoMau(pink);
		base.Draw(mode, amb, spec, shi);
		glPopMatrix();

		for (int i = 0; i <= 5; i++) {
			glPushMatrix();
			glRotatef(60 * i, 0, 1, 0);
			glTranslatef(1.4, 3.6, 0);
			khoiBanNguyet.KhoiTaoMau(blue);
			khoiBanNguyet.Draw(mode, amb, spec, shi);
			glPopMatrix();
		}
		glPopMatrix();
	}
	else {
		glDisable(GL_LIGHTING);
		glPushMatrix();
		glTranslatef(7.2, 0, 0);
		glRotatef(xoay, 0, 1, 0);
		glTranslatef(-7.2, 0, 0);
		glTranslatef(3.6, 6.6, 0);
		hinhTru.DrawKhung();
		glPopMatrix();
		glPushMatrix();
		glTranslatef(7.2, 5.8, 0);
		glRotatef(xoay, 0, 1, 0);
		khoiHopXoay.DrawKhung();
		glPopMatrix();
		glPushMatrix();
		glTranslatef(7.2, 2.2 + 1.5 * 0.5 + 3.6, 0);
		glRotatef(xoay, 0, 1, 0);
		hinhHop.DrawKhung();
		glPopMatrix();
		glPushMatrix();
		glMultMatrixf(modelview);
		glPushMatrix();
		glTranslatef(0, 3.6 * 0.5, 0);
		base.DrawKhung();
		glPopMatrix();

		for (int i = 0; i <= 5; i++) {
			glPushMatrix();
			glRotatef(60 * i, 0, 1, 0);
			glTranslatef(1.4, 3.6, 0);
			khoiBanNguyet.DrawKhung();
			glPopMatrix();
		}
		glPopMatrix();
		glEnable(GL_LIGHTING);
	}
}

void create() {
	khoiBanNguyet.KhoiBanNguyet(28, 2.2, 3.8, sqrt(3.0) * 5.8 / 2, 5.6, 12.5, 60, 65, 5.8);
	khoiBanNguyet.TinhGocQuay();
	khoiHopXoay.KhoiHopXoay(28, 1.5, 2.6, 4.8, 5.6, 72);
	khoiHopXoay.TinhGocQuay();
	hinhHop.HinhHop(7.2, 1.5, 1.4);
	hinhHop.TinhGocQuay();
	hinhTru.HinhTru(28, 6.0, 0.7);
	hinhTru.TinhGocQuay();
	base.HinhHop(31.5, 3.6, 31.5);
	base.TinhGocQuay();
}

void display() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	view2D ? glOrtho(-33.0, 33.0, -33.0 / 2, 33.0 / 2, -100.0, 100.0) : gluPerspective(30.0, (float)scrWidth / (float)scrHeight, 1.0, 10000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	if(!view2D) {
		(camDis == 0) ? gluLookAt(camDis * sinf(camGocQuay * PI / 180), camChieuRong, camDis * cosf(camGocQuay * PI / 180), 0, 4, 0, sinf(camGocQuay * PI / 180), 0, cosf(camGocQuay * PI / 180))
			: gluLookAt(camDis * sinf(camGocQuay * PI / 180), camChieuRong, camDis * cosf(camGocQuay * PI / 180), 0, 2, 0, 0, 1, 0);
	} else {
		gluLookAt(0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	}

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, scrWidth, scrHeight);
	GLfloat mat_diff[4] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_spec[4] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_amb[4] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lightPosition_0[4] = { -9.0, 15.0, 11.0, 0.0 };
	GLfloat lightPosition_1[4] = { 9.0, 15.0, -11.0, 0.0 };

	glLightfv(GL_LIGHT0, GL_DIFFUSE, mat_diff);
	glLightfv(GL_LIGHT0, GL_AMBIENT, mat_amb);
	glLightfv(GL_LIGHT0, GL_SPECULAR, mat_spec);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition_0);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, mat_diff);
	glLightfv(GL_LIGHT1, GL_AMBIENT, mat_amb);
	glLightfv(GL_LIGHT1, GL_SPECULAR, mat_spec);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPosition_1);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	drawNen(0.6f);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glPushMatrix();
	glTranslatef(0, -2, 0);
	drawAllMesh();
	glPopMatrix();
	glFlush();
	glutSwapBuffers();
}

void myInit() {
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glFrontFace(GL_CCW);
	glEnable(GL_NORMALIZE);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
}

void printUsage() {
	cout << "1		: Quay nguoc chieu kim dong ho" << endl;
	cout << "2		: Quay cung chieu kim dong ho" << endl;
	cout << "3		: Bat/tat che do to mau nang cao" << endl;
	cout << "V, v		: Chuyen doi giua 2 che do nhin" << endl;
	cout << "W, w		: Chuyen doi qua lai giua che do khung day va to mau" << endl;
	cout << "+		: Tang khoang cach camera" << endl;
	cout << "-		: Giam khoang cach camera" << endl;
	cout << "up arrow	: Tang chieu cao camera" << endl;
	cout << "down arrow	: Giam chieu cao camera" << endl;
	cout << "<-		: Quay camera theo chieu kim dong ho" << endl;
	cout << "->		: Quay camera nguoc chieu kim dong ho" << endl;
}

int main(int argc, char* argv[]) {
	printUsage();
	glutInit(&argc, (char**)argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(scrWidth, scrHeight);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Che Quang Huy - 1812340");
	create();
	myInit();
	glutDisplayFunc(display);
	glutSpecialFunc(controlKeyboard);
	glutKeyboardFunc(control);
	glutMainLoop();
	return 0;
}