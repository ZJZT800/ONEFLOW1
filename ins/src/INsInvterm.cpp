/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

//#include "UINsCorrectSpeed.h"
#include "INsInvterm.h"
#include "INsVisterm.h"
#include "UINsCom.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "Com.h"
#include "INsCom.h"
#include "INsIDX.h"
#include "HXMath.h"
#include "Ctrl.h"

BeginNameSpace( ONEFLOW )

INsInv iinv;



INsInv::INsInv()
{
    ;
}

INsInv::~INsInv()
{
    ;
}

void INsInv::Init()
{
    int nEqu = inscom.nEqu;
    prim.resize( nEqu );
    prim1.resize( nEqu );
    prim2.resize( nEqu );

    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

    dq.resize( nEqu );

    flux.resize( nEqu );
    flux1.resize( nEqu );
    flux2.resize( nEqu );

}

INsInvterm::INsInvterm()
{
    ;
}

INsInvterm::~INsInvterm()
{
    ;
}

void INsInvterm::Solve()
{
}

void INsInvterm::CmpINsinvTerm()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	Real v2l = ONEFLOW::SQR(iinv.ul, iinv.vl, iinv.wl);
	Real v2r = ONEFLOW::SQR(iinv.ur, iinv.vr, iinv.wr);

	Real vnl = gcom.xfn * iinv.ul + gcom.yfn * iinv.vl + gcom.zfn * iinv.wl - gcom.vfn;       // V * n
	Real vnr = gcom.xfn * iinv.ur + gcom.yfn * iinv.vr + gcom.zfn * iinv.wr - gcom.vfn;

	Real rvnl = iinv.rl * vnl;   //�� * V * n
	Real rvnr = iinv.rr * vnr;


	iinv.rf[ug.fId] = (iinv.rl + iinv.rr) * half;    //��ʼ�����ϵ�ֵ��u��v��w ��
	iinv.uf[ug.fId] = (iinv.ul + iinv.ur ) * half;
	iinv.vf[ug.fId] = (iinv.vl + iinv.vr ) * half;
	iinv.wf[ug.fId] = (iinv.wl + iinv.wr ) * half;

	
	iinv.vnflow[ug.fId] = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId];  //��ʼ������ V*n

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow[ug.fId] * gcom.farea; //��ʼ�����ϵ�����ͨ��

     Real clr = MAX(0, iinv.fq[ug.fId]);  //�ӽ�����൥Ԫ�����Ҳ൥Ԫ�ĳ�ʼ��������
     Real crl = clr - iinv.fq[ug.fId];   //�ӽ����Ҳ൥Ԫ������൥Ԫ�ĳ�ʼ��������


	iinv.aii1[ug.lc] = crl;   //����������Ԫ������
	iinv.aii2[ug.rc] = clr;   //���������ҵ�Ԫ������
	
    iinv.ai1[ug.lc]+= crl;   //���뵥Ԫ������
	iinv.ai2[ug.rc]+=clr;   //������Ԫ������
	
}

void INsInvterm::CmpINsFaceflux()
{
	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
	Real dy2 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
	Real dz2 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

	Real de1 = DIST(dx1, dy1, dz1);
	Real de2 = DIST(dx2, dy2, dz2);
	Real de = 1.0 / (de1 + de2);

	iinv.f1[ug.fId] = de2 * de;  //��ԪȨ��
    iinv.f2[ug.fId] = de1 * de;  //�ҵ�ԪȨ��
 
	Real Vau = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spu[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spu[ug.rc]);  //��Vj/a��
	Real Vav = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spv[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spv[ug.rc]);
	Real Vaw = iinv.f1[ug.fId] * (gcom.cvol1 / iinv.spw[ug.lc]) + iinv.f2[ug.fId] * (gcom.cvol2 / iinv.spw[ug.rc]);

    iinv.dist[ug.fId] = gcom.xfn * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);
	Real Pd1 = visQ.dqdx1[IIDX::IIP] * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + visQ.dqdy1[IIDX::IIP] * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + visQ.dqdz1[IIDX::IIP] * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);  //ѹ���ݶ���
	Real Pd2 = visQ.dqdx2[IIDX::IIP] * ((*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId]) + visQ.dqdy2[IIDX::IIP] * ((*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId]) + visQ.dqdz2[IIDX::IIP] * ((*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId]);
	Real Pd = Pd1 + Pd2;

	
	iinv.uf[ug.fId] = (iinv.f1[ug.fId] *iinv.ul + iinv.f2[ug.fId] *iinv.ur) + (iinv.Vau*gcom.xfn / iinv.dist[ug.fId])*( Pd - (iinv.pr - iinv.pl));  //��һʱ�̵Ľ���Ԥ���ٶ�
	iinv.vf[ug.fId] = (iinv.f1[ug.fId] *iinv.vl + iinv.f2[ug.fId] *iinv.vr) + (iinv.Vav*gcom.yfn / iinv.dist[ug.fId])*(Pd - (iinv.pr - iinv.pl));
	iinv.wf[ug.fId] = (iinv.f1[ug.fId] *iinv.wl + iinv.f2[ug.fId] *iinv.wr) + (iinv.Vaw*gcom.zfn / iinv.dist[ug.fId])*(Pd - (iinv.pr - iinv.pl));

	iinv.vnflow[ug.fId] = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId];
	

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow[ug.fId] * gcom.farea;  //��һʱ�̽���Ԥ��ͨ��
}

void INsInvterm::CmpINsFaceCorrectPresscoef()
{

	iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ( gcom.cvol1/((1+1)*iinv.spu[ug.lc] - iinv.sj[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 2)*iinv.spu[ug.rc] - iinv.sj[ug.rc]));  // (Vp/dv)j�����������ٶ�������
	iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * (gcom.cvol1 / ((1+1)*iinv.spv[ug.lc] - iinv.sj[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 2)*iinv.spv[ug.rc] - iinv.sj[ug.rc]));
	iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * (gcom.cvol1 / ((1+1)*iinv.spw[ug.lc] - iinv.sj[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 2)*iinv.spw[ug.rc] - iinv.sj[ug.rc]));
	

	iinv.aju[ug.fId] = iinv.rm * iinv.Vdvu[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId]; //ѹ�����������У������з���ϵ��
	iinv.ajv[ug.fId] = iinv.rm * iinv.Vdvv[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId];
	iinv.ajw[ug.fId] = iinv.rm * iinv.Vdvw[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId];

	iinv.ajp[ug.fId] = iinv.aju[ug.fId] * gcom.xfn + iinv.ajv[ug.fId] * gcom.yfn + iinv.ajw[ug.fId] * gcom.zfn;
}


  














EndNameSpace