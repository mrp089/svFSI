// NOTE: This file is automatically included from tens4d.h
// Users should not include this file manually!

#include "matrix.h"

inline tens4dss::tens4dss(const double g)
{
    for (int i = 0; i < NNZ; i++)
        d[i] = g;
}

inline tens4dss::tens4dss(double m[6][6])
{
    d[0] = m[0][0]; d[ 6] = m[0][1]; d[12] = m[0][2]; d[18] = m[0][3]; d[24] = m[0][4]; d[30] = m[0][5];
    d[1] = m[1][0]; d[ 7] = m[1][1]; d[13] = m[1][2]; d[19] = m[1][3]; d[25] = m[1][4]; d[31] = m[1][5];
    d[2] = m[2][0]; d[ 8] = m[2][1]; d[14] = m[2][2]; d[20] = m[2][3]; d[26] = m[2][4]; d[32] = m[2][5];
    d[3] = m[3][0]; d[ 9] = m[3][1]; d[15] = m[3][2]; d[21] = m[3][3]; d[27] = m[3][4]; d[33] = m[3][5];
    d[4] = m[4][0]; d[10] = m[4][1]; d[16] = m[4][2]; d[22] = m[4][3]; d[28] = m[4][4]; d[34] = m[4][5];
    d[5] = m[5][0]; d[11] = m[5][1]; d[17] = m[5][2]; d[23] = m[5][3]; d[29] = m[5][4]; d[35] = m[5][5];
}

inline tens4dss::tens4dss(tens4ds& t)
{
	d[0] = t.d[ 0]; d[ 6] = t.d[ 1]; d[12] = t.d[ 3]; d[18] = t.d[ 6]; d[24] = t.d[10]; d[30] = t.d[15];
    d[1] = t.d[ 1]; d[ 7] = t.d[ 2]; d[13] = t.d[ 4]; d[19] = t.d[ 7]; d[25] = t.d[11]; d[31] = t.d[16];
    d[2] = t.d[ 3]; d[ 8] = t.d[ 4]; d[14] = t.d[ 5]; d[20] = t.d[ 8]; d[26] = t.d[12]; d[32] = t.d[17];
    d[3] = t.d[ 6]; d[ 9] = t.d[ 7]; d[15] = t.d[ 8]; d[21] = t.d[ 9]; d[27] = t.d[13]; d[33] = t.d[18];
    d[4] = t.d[10]; d[10] = t.d[11]; d[16] = t.d[12]; d[22] = t.d[13]; d[28] = t.d[14]; d[34] = t.d[19];
    d[5] = t.d[15]; d[11] = t.d[16]; d[17] = t.d[17]; d[23] = t.d[18]; d[29] = t.d[19]; d[35] = t.d[20];
}

inline double& tens4dss::operator () (int i, int j, int k, int l)
{
    const int m[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
    tens4dss& T = (*this);
    return T(m[i][j], m[k][l]);
}

inline double tens4dss::operator () (int i, int j, int k, int l) const
{
    const int m[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
    const tens4dss& T = (*this);
    return T(m[i][j], m[k][l]);
}

inline double& tens4dss::operator () (int i, int j)
{
    const int m[6] = {0, 6, 12, 18, 24, 30};
    return d[m[j]+i];
}

inline double tens4dss::operator () (int i, int j) const
{
    const int m[6] = {0, 6, 12, 18, 24, 30};
    return d[m[j]+i];
}

// operator +
inline tens4dss tens4dss::operator + (const tens4dss& t) const
{
    tens4dss s;
    for (int i=0; i<NNZ; i++)
        s.d[i] = d[i] + t.d[i];
    
    return s;
}

// operator -
inline tens4dss tens4dss::operator - (const tens4dss& t) const
{
    tens4dss s;
    for (int i=0; i<NNZ; i++)
        s.d[i] = d[i] - t.d[i];

    return s;
}

// operator *
inline tens4dss tens4dss::operator * (double g) const
{
    tens4dss s;
    for (int i=0; i<NNZ; i++)
        s.d[i] = g*d[i];
    
    return s;
}

// operator /
inline tens4dss tens4dss::operator / (double g) const
{
    tens4dss s;
    for (int i=0; i<NNZ; i++)
        s.d[i] = d[i]/g;
    
    return s;
}

// assignment operator +=
inline tens4dss& tens4dss::operator += (const tens4dss& t)
{
    for (int i=0; i<NNZ; i++)
        d[i] += t.d[i];
    
    return (*this);
}

// assignment operator -=
inline tens4dss& tens4dss::operator -= (const tens4dss& t)
{
    for (int i=0; i<NNZ; i++)
        d[i] -= t.d[i];
    
    return (*this);
}

// assignment operator *=
inline tens4dss& tens4dss::operator *= (double g)
{
    for (int i=0; i<NNZ; i++)
        d[i] *= g;
    
    return (*this);
}

// assignment operator /=
inline tens4dss& tens4dss::operator /= (double g)
{
    for (int i=0; i<NNZ; i++)
        d[i] /= g;
    
    return (*this);
}

// unary operator -
inline tens4dss tens4dss::operator - () const
{
    tens4dss s;
    for (int i = 0; i < NNZ; i++)
        s.d[i] = -d[i];

    return s;
}

// intialize to zero
inline void tens4dss::zero()
{
    for (int i = 0; i < NNZ; i++)
        d[i] = 0;
}

// extract 6x6 matrix
inline void tens4dss::extract(double D[6][6])
{
    D[0][0] = d[0]; D[0][1] = d[ 6]; D[0][2] = d[12]; D[0][3] = d[18]; D[0][4] = d[24]; D[0][5] = d[30];
    D[1][0] = d[1]; D[1][1] = d[ 7]; D[1][2] = d[13]; D[1][3] = d[19]; D[1][4] = d[25]; D[1][5] = d[31];
    D[2][0] = d[2]; D[2][1] = d[ 8]; D[2][2] = d[14]; D[2][3] = d[20]; D[2][4] = d[26]; D[2][5] = d[32];
    D[3][0] = d[3]; D[3][1] = d[ 9]; D[3][2] = d[15]; D[3][3] = d[21]; D[3][4] = d[27]; D[3][5] = d[33];
    D[4][0] = d[4]; D[4][1] = d[10]; D[4][2] = d[16]; D[4][3] = d[22]; D[4][4] = d[28]; D[4][5] = d[34];
    D[5][0] = d[5]; D[5][1] = d[11]; D[5][2] = d[17]; D[5][3] = d[23]; D[5][4] = d[29]; D[5][5] = d[35];
}

//-----------------------------------------------------------------------------
// (a dyad1ss b)_ijkl = a_ij b_kl
inline tens4dss dyad1ss(const mat3ds& a, const mat3ds& b)
{
    tens4dss c;
    
    c.d[ 0] = a.xx()*b.xx();
    c.d[ 1] = a.yy()*b.xx();
    c.d[ 2] = a.zz()*b.xx();
    c.d[ 3] = a.xy()*b.xx();
    c.d[ 4] = a.yz()*b.xx();
    c.d[ 5] = a.xz()*b.xx();
    
    c.d[ 6] = a.xx()*b.yy();
    c.d[ 7] = a.yy()*b.yy();
    c.d[ 8] = a.zz()*b.yy();
    c.d[ 9] = a.xy()*b.yy();
    c.d[10] = a.yz()*b.yy();
    c.d[11] = a.xz()*b.yy();
    
    c.d[12] = a.xx()*b.zz();
    c.d[13] = a.yy()*b.zz();
    c.d[14] = a.zz()*b.zz();
    c.d[15] = a.xy()*b.zz();
    c.d[16] = a.yz()*b.zz();
    c.d[17] = a.xz()*b.zz();
    
    c.d[18] = a.xx()*b.xy();
    c.d[19] = a.yy()*b.xy();
    c.d[20] = a.zz()*b.xy();
    c.d[21] = a.xy()*b.xy();
    c.d[22] = a.yz()*b.xy();
    c.d[23] = a.xz()*b.xy();
    
    c.d[24] = a.xx()*b.yz();
    c.d[25] = a.yy()*b.yz();
    c.d[26] = a.zz()*b.yz();
    c.d[27] = a.xy()*b.yz();
    c.d[28] = a.yz()*b.yz();
    c.d[29] = a.xz()*b.yz();
    
    c.d[30] = a.xx()*b.xz();
    c.d[31] = a.yy()*b.xz();
    c.d[32] = a.zz()*b.xz();
    c.d[33] = a.xy()*b.xz();
    c.d[34] = a.yz()*b.xz();
    c.d[35] = a.xz()*b.xz();
    
    return c;
}

//-----------------------------------------------------------------------------
// (a ddotss b)_ijkl = a_ijmn b_mnkl
inline tens4dss ddotss(const tens4dss& a, const tens4dss& b)
{
	tens4dss c;
	
	// compute c in matrix notation
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			c(i,j) = (a(i,0)*b(0,j) + a(i,1)*b(1,j) + a(i,2)*b(2,j))
				   + (a(i,3)*b(3,j) + a(i,4)*b(4,j) + a(i,5)*b(5,j))*2.0;
		}
	}
    
    return c;
}
