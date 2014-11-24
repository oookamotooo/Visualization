#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "get_jacob.h"
#include "Vector3.h"
#include <vector>
#include "define.h"

//global define

int r_number0;
int r_number1;
int r_number2;


double arc;
double L = 0;
double add = 0.5;

double r_delta_max;
double band_delta_min;

//simulation parameter
double Delta = 0.1;
int split = 100;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���s�O��step�̒l�m�F!!!!!!!!!!!!!!!!!!!!!

double step = -0.1;
double runtime_out = 150;
double len = 0.01;

//�v���g�^�C�v�錾(����Ƃ��s�v�Ȃ��)
/*
int deltalength_matching(double band_delta_min, int r_number2,  double Band1[r_number1][3], double Band2[r_number2][3]);
void read_data();
*/

//BVP�֐�
Vector3d BVP(Vector3d B0[band_size], Vector3d B1[band_size], Vector3d leaf, int r, int cpnum, double Delta, Vector3d critical[CP_NUM]){
    int i;
    int informal_r;
    int runtime = 0;//�����Q�N�b�^������J��Ԃ�����runtime_out
    int taw;
    double tawsp=0;
    int discover = 0;
    double r_distance;
    double br_distance;
    double difference;
    double naiseki;
    double yobidelta;
    double rtoe_br_length;
    double mukinaiseki;
	double isFar_1, isFar_r;
	double var_size;
	Vector3d tmp;
	bool noFar;

	Vector3d br, informal_br, e_br, hosei, Vec01;
	Vector3d solution;
	/*
	double br[3];
    double informal_br[3];
    double e_br[3];
	double hosei[3];
	double Vec01[3];
    */
	/*
    FILE *fp;
    if (NULL == (fp = fopen("./debag_solution.txt", "w"))) {
        fprintf(stderr, "cannot open output file\n");
        exit(-1);
    }*/
    
	//������
	br = Vector3d();	informal_br = Vector3d();	e_br = Vector3d();
    yobidelta = 10;	noFar = true;
    
	printf("First noFar = %d\n", noFar);
    //Band1[0]���珇�ɃV���[�e�B���O����
    for(informal_r = 0; informal_r < r_number1-1; informal_r++){
        for(taw=0;taw<=split;taw++){
            
			//���_����������1
            discover = 0;
            
            tawsp = ((double)taw/(double)split);
            /* if(informal_r == r_number1-1){
             br[0] = (double)(Band1[informal_r][0]*(1.0-tawsp) + Band1[0][0]*tawsp);
             br[1] = (double)(Band1[informal_r][1]*(1.0-tawsp) + Band1[0][1]*tawsp);
             br[2] = (double)(Band1[informal_r][2]*(1.0-tawsp) + Band1[0][2]*tawsp);
             //printf("br if !\n");
             }else{
            br[0] = (double)(Band1[informal_r][0]*(1.0-tawsp) + Band1[informal_r+1][0]*tawsp);
            br[1] = (double)(Band1[informal_r][1]*(1.0-tawsp) + Band1[informal_r+1][1]*tawsp);
            br[2] = (double)(Band1[informal_r][2]*(1.0-tawsp) + Band1[informal_r+1][2]*tawsp);
            */
			br = B1[informal_r]*(1.0-tawsp) + B1[informal_r+1]*tawsp;

            runtime = 0;

            while(runtime < runtime_out){
               
				//�����Q�N�b�^����
				/*
				tmp.x = runge_kutta(br.x, br.y, br.z, 1, step);
				tmp.y = runge_kutta(br.x, br.y, br.z, 2, step);
				tmp.z = runge_kutta(br.x, br.y, br.z, 3, step);
				*/
				var_size = Vector3d().dist(tmp);

				br.x = br.x + (tmp.x/var_size)*len;
				br.y = br.y + (tmp.y/var_size)*len;
				br.z = br.z + (tmp.z/var_size)*len;

				/*
                informal_br[0] = (br[0] - Band1[r][0]);
                informal_br[1] = (br[1] - Band1[r][1]);
                informal_br[2] = (br[2] - Band1[r][2]);
                */
                informal_br = br - B1[r];

                naiseki = 0;
                //naiseki = (double)(informal_br.x*leaf.x + informal_br.y*leaf[1] + informal_br[2]*leaf[2]);
				naiseki = informal_br.dot(leaf);
                
                difference = 0;
                //difference = fabs((double)(sqrt(pow((br[0]-B1[r][0]), 2) + pow((br[1]-B1[r][1]), 2) + pow((br[2]-B1[r][2]), 2))));
				difference = fabs(br.dist(B1[r]));
				
				/*
                r_distance = 0;
                br_distance = 0;
                r_distance = fabs((double)(sqrt(pow(Band1[r][0], 2) + pow(Band1[r][1], 2) + pow(Band1[r][2], 2))));
                br_distance = fabs((double)(sqrt(pow(br[0], 2) + pow(br[1], 2) + pow(br[2], 2))));
                //difference-Delta�����łȂ���Ύ��Ȃ��悤�ɂ��Ă݂�H
                */
				
				//!!!!!!�����o�[�W�����ɂ��邽�ߍ폜
				Vec01 = B1[r] - B0[r];
				mukinaiseki = informal_br.dot(Vec01);


				//!!!!!!!!!�����ꍇ�͕K���������ٓ_���牓����΂������Ă킯�ł��Ȃ�!!!!!
				/*isFar_1 = B1[r].dist(critical[cpnum]);
                isFar_r = br.dist(critical[cpnum]);*/

                if(fabs(naiseki) < 0.1){//leaf���`����@���x�N�g���ƒ��p�ł��邩�ǂ���
					//printf("naiseki\n");
                    if(/*(isFar_r - isFar_1) > 0*/mukinaiseki > 0){//Vec01(Band1-B0)�Ɛ����p�x�@�s�p�ł���΂悵 !!!!����cp to B1 ��cp to br�̋�����r����΂����񂶂�Ȃ���
                        //printf("isFar\n");
						noFar = false;
						if(r_delta_max < difference){
                            r_delta_max = difference;
                        }
                        if(fabs(difference-Delta)<yobidelta){//���z�Ƃ��鋗�����ɋ߂����ǂ���
                            yobidelta = fabs(difference-Delta);
                            /*
                            e_br[0] = br[0];
                            e_br[1] = br[1];
                            e_br[2] = br[2];
                            */
							e_br = br;
                        }
                        if(fabs(difference-Delta)<0.05){
							//printf("Find!\n");
                            //�K���������_����̋�����������ΐVbr�Ƃ��Đ������Ƃ͌����Ȃ��̂ŉ��ǂ̗]�n����!!
                            /*
							solution[0] = br[0];
                            solution[1] = br[1];
                            solution[2] = br[2];
                            */
							solution = br;
                            discover = 1;
                            //printf("discover new Point!\n");
                            break;
                        }
                    }
                }
                runtime++;
            }
            if(discover == 1)break;
        }
        if(discover == 1)break;
    }
    
	printf("noFar = %d\n", noFar);
	if(noFar){
		solution = B1[r];
		printf("The point far from B1 is not found.\n");
	}else if(discover == 0){
        
        //rtoe_br_length = sqrt(pow((e_br[0]-B1[r][0]), 2) + pow((e_br[1]-B1[r][1]), 2) + pow((e_br[2]-B1[r][2]), 2));
		rtoe_br_length = e_br.dist(B1[r]);
		
		/*
		hosei[0] = (e_br[0]-B1[r][0])*(Delta/rtoe_br_length);
        hosei[1] = (e_br[1]-B1[r][1])*(Delta/rtoe_br_length);
        hosei[2] = (e_br[2]-B1[r][2])*(Delta/rtoe_br_length);
        */
		hosei = (e_br - B1[r])*(Delta/rtoe_br_length);

		/*
        solution[0] = (B1[r][0] + hosei[0]);
        solution[1] = (B1[r][1] + hosei[1]);
        solution[2] = (B1[r][2] + hosei[2]);
		*/
		solution = (B1[r] + hosei);
		printf("Best point not found\n");
    }
    //�V���[�e�B���O�œ����ꂽ�e�_�Ɋւ���r��br�x�N�g����leaf�����s���Ă邩���ׂ�
    //���s���Ă��Ƃ��ċ��������ł��肩��r��艓���ʒu�ɂ��邩�ǂ�������
    //�����𖞂����Ă����_��solution�Ƃ��đ��
    
    //fclose(fp);
    return solution;
}

//Accuracy�֐�
int Accuracy_long(Vector3d AB0[band_size], Vector3d AB1[band_size], Vector3d AB2[band_size], int cpnum, Vector3d critical[CP_NUM]){
    int i, s, l;
    double mesh_length;
    Vector3d leaf;
    
    for(i=0;i<r_number2-1;i++){
        
        //mesh_length = (double)(sqrt(pow((AB2[i][0]-AB2[i+1][0]), 2) + pow((AB2[i][1]-AB2[i+1][1]), 2) + pow((AB2[i][2]-AB2[i+1][2]), 2)));
		mesh_length = AB2[i].dist(AB2[i+1]);

        if(mesh_length > 2.0){
            r_number1++;
            r_number2++;
            
            for(s=r_number1-1;s>i;s--){
				/*
                AB1[s][0] = AB1[s-1][0];
                AB1[s][1] = AB1[s-1][1];
                AB1[s][2] = AB1[s-1][2];
                
                AB2[s][0] = AB2[s-1][0];
                AB2[s][1] = AB2[s-1][1];
                AB2[s][2] = AB2[s-1][2];
                */
				AB1[s] = AB1[s-1];
				AB2[s] = AB2[s-1];
            }
            /*
            AB1[i+1][0] = (AB1[i][0] + AB1[i+2][0])/2.0;
            AB1[i+1][1] = (AB1[i][1] + AB1[i+2][1])/2.0;
            AB1[i+1][2] = (AB1[i][2] + AB1[i+2][2])/2.0;
            */
			AB1[i+1] = (AB1[i] + AB1[i+2])/2.0;

			/*
            leaf[0] = AB1[i+2][0] - AB1[i][0];
            leaf[1] = AB1[i+2][1] - AB1[i][1];
            leaf[2] = AB1[i+2][2] - AB1[i][2];
            */
			leaf = AB1[i+2] - AB1[i];

			/*
            solution[0] = 0;
            solution[1] = 0;
            solution[2] = 0;
			*/
            
            AB2[i+1] = BVP(AB0, AB1, leaf, i+1, cpnum, Delta, critical);
            
            if(band_delta_min > r_delta_max){
                band_delta_min = r_delta_max;
            }
            /*
            AB2[i+1][0] = solution[0];
            AB2[i+1][1] = solution[1];
            AB2[i+1][2] = solution[2];
            */

            /*for(l=0;l<r_number2;l++){
             printf("Band2[%d] = %f %f %f\n",l, Band2[l][0], Band2[l][1], Band2[l][2] );
             }*/
            
        }
    }
    return 1;
}

int Accuracy_short(Vector3d AB2[band_size]){
    int i, s;
    double mesh_length;
    
    for(i=0;i<r_number2-1;i++){
        
        //mesh_length = (double)(sqrt(pow((AB2[i][0]-AB2[i+1][0]), 2) + pow((AB2[i][1]-AB2[i+1][1]), 2) + pow((AB2[i][2]-AB2[i+1][2]), 2)));
        mesh_length = AB2[i].dist(AB2[i+1]);

        if(mesh_length < 0.5){
            printf("---------------mesh_length is short!-------------\n");
            printf("mesh_length is %f !\n", mesh_length);
            
            for(s=i;s<r_number2-1;s++){
				/*
                AB2[s][0] = AB2[s+1][0];
                AB2[s][1] = AB2[s+1][1];
                AB2[s][2] = AB2[s+1][2];
				}*/
				AB2[s] = AB2[s+1];
            }
            /*
            AB2[r_number2-1][0] = 0;
            AB2[r_number2-1][1] = 0;
            AB2[r_number2-1][2] = 0;
            */
			AB2[r_number2-1] = Vector3d();

            r_number2--;
        }
    }
    return 1;
}

//main�֐�
void Geodesic(Vector3d critical[CP_NUM], Vector3d critical_round[CP_NUM][round_num], int cpnum){
    //int r_number1 = 10;
    //double arc, L=5;
    int i, r;

    int Acc_l, Acc_s, delmatch = 0;
    Vector3d Band0[band_size];
    Vector3d Band1[band_size];
    Vector3d Band2[band_size];

    Vector3d leaf;
    
    int sub_number0, sub_number1;
    
    int Admissible = 0;

	FILE *fp;
	//char name = (char)n;
	char filepath[256];
	
	sprintf(filepath, "cp%d_flow_noFar2.txt", cpnum);
	fp = fopen(filepath, "w");
	
	/*
    FILE *fp;
    if (NULL == (fp = fopen("./lorenz_manifold90_91_delta.txt", "w"))) {
        fprintf(stderr, "cannot open output file\n");
        exit(-1);
    }
    FILE *fp2;
    if (NULL == (fp2 = fopen("./savedata91_delta.txt", "w"))) {
        fprintf(stderr, "cannot open output file\n");
        exit(-1);
    }*/
    
	//fprintf(fp, "\n");
    
	r_number0 = round_num;	r_number1 = round_num; r_number2 = round_num;

	//������Band�Ɋi�[
	for(int i = 0; i<r_number1; i++){
		Band0[i] = critical[cpnum];
		Band1[i] = critical[cpnum] + critical_round[cpnum][i];
		fprintf(fp, "%lf\t%lf\t%lf\n", Band1[i].x, Band1[i].y, Band1[i].z);
	}

    //arc = sqrt(pow(Band1[0], 2) + pow(Band1[0], 2) + pow(Band1[0][2], 2));
	arc = Vector3d().dist(Band1[0]);

	printf("First arc length is %f\n", arc);
    printf("ababa");
    
    L = arc + add;
    printf("L = %f\n",L);
    
    while(arc < L){ //�ʒ�L�ɂȂ�܂ő�����
        //band_delta_min��������
        r_delta_max = 0;
        band_delta_min = 1.0;
		Acc_l = 0;
        Admissible = 0;

		printf("r_number1 = %d\n", r_number1);
        for(r=0;r<r_number1;r++){
            //�����Ă���Band1�̊e�_r�ɂ��ėt�w�\����`���V���[�e�B���O�@��Band2��r��`�����b�V���̒ǉ��폜��Band1��Band0�ɁABand2��Band1�Ɂ@���J��Ԃ��B
            
			//round_cp�Ɋւ��Č�����if���ł�0�A�h���X�ꍇ�킯�K�v�Ȃ��H���� 
			/*
			leaf[0] = Band1[r+1][0] - Band1[r-1][0];
			leaf[1] = Band1[r+1][1] - Band1[r-1][1];
			leaf[2] = Band1[r+1][2] - Band1[r-1][2];
            */
			if(r == 0){
				leaf = Band1[r+1] - Band1[r_number1-1];
			}else if(r == r_number1-1){
				leaf = Band1[0] - Band1[r-1];
			}else{
				leaf = Band1[r+1] - Band1[r-1];
			}

            //�V�_(Cb, Band2)��`
            /*
			solution[0] = 0;
			solution[1] = 0;
			solution[2] = 0;
			*/
			
			//solution = Vector3d();
                
			Band2[r] = BVP(Band0, Band1, leaf, r, cpnum, Delta, critical);
                
			//BVP�ɂ���č����Ă��邒�ɂ��Ă̑������ʂ���Aleaf��̂����Ƃ������_�����߂�ꂽ�B���ꂪ�����̒��ōł��Z���ꍇ�A�����ێ�
            //�����͌�X�ł���
			/*
			if(band_delta_min > r_delta_max){
				band_delta_min = r_delta_max;
				printf("band delta min = %f\n", band_delta_min);
			}*/

            /*    
			Band2[r][0] = solution[0];
			Band2[r][1] = solution[1];
			Band2[r][2] = solution[2];
			*/

            printf("r = %d", r);
			printf(":%lf %lf %lf\n", Band2[r].x, Band2[r].y, Band2[r].z);
	}
        
		//���b�V���Ԃ��������Ȃ���check
        Acc_l = Accuracy_long(Band0, Band1, Band2, cpnum, critical);
        
		//delmatch���Ƃ�
		/*
        delmatch = 0;
        if(Delta > band_delta_min){
            delmatch = deltalength_matching(band_delta_min,r_number2,Band1,Band2);//�~�̐L�ї�Delta��(�K�v�ł����)�␳
        }*/
        
		//���b�V���Ԃ��Z�����Ȃ���check
        //Acc_s = Accuracy_short(Band2);
        
        //��̔���֐��Ŗ��Ȃ����Admissible��1��
        if(/*Acc_s == 1 && */Acc_l == 1){
            Admissible = 1;
        }
        
		//�ȉ�Band�̍X�V��fprintf�ɂ�鏑���o��
        if(Admissible){
            if(delmatch == 1){
                printf("delmatch == 1! Arclength = %f ��", arc);
                arc = arc + band_delta_min;
                printf("%f\n", arc);
            }else{
                arc = arc + Delta;
                printf("Arclength = %f\n", arc);
            }
            
            for(r=0;r<r_number1;r++){
                //Band0[r][0] = Band1[r][0]; Band0[r][1] = Band1[r][1]; Band0[r][2] = Band1[r][2];
				Band0[r] = Band1[r];
            }
            for(r=0;r<r_number1;r++){
                //Band1[r][0] = 0; Band1[r][1] = 0; Band1[r][2] = 0;
				Band1[r] = Vector3d();
            }
            for(r=0;r<r_number2;r++){
                //Band1[r][0] = Band2[r][0]; Band1[r][1] = Band2[r][1]; Band1[r][2] = Band2[r][2];
				Band1[r] = Band2[r];
                
				printf("r%d\n", r);
                //Accuracy�ŗv�f���ɍ����o���ꍇ�ǂ����邩�H���Ƃ��Ă�Band������������ʂɊm�ۂ��Ƃ��āA�Q�Ƃ͊��S��r_number1�ɗ������
                //fprintf(fp, "%lf\t%lf\t%lf\n", Band0[r][0], Band0[r][1], Band0[r][2]);
                fprintf(fp, "%lf\t%lf\t%lf\n", Band1[r].x, Band1[r].y, Band1[r].z);
                //printf("Cr=%d : %f %f %f\n", r, Band1[r][0], Band1[r][1], Band1[r][2]);
            }
			/*
            for(r=1;r<r_number2-1;r++){
                fprintf(fp, "%lf\t%lf\t%lf\n", -Band1[r_number2-1-r]);
                //printf("Cr=%d : %f %f %f\n", r_number2-1+r, -Band1[r_number2-1-r][0], -Band1[r_number2-1-r][1], Band1[r_number2-1-r][2]);
            }*/


            sub_number0 = r_number1;
            r_number1 = r_number2;
            sub_number1 = r_number2;
        }else{
				printf("!!!!!!!!!!!!!!!!!!Admissible 0!!!!!!!!!!!!!!!!!\n");
		}
        //fprintf(fp, "\n");
        //adjust Delta Algorithm
    }
    
    //�����܂�while
    /*
    fprintf(fp2, "%d\n", sub_number0);
    for(r=0;r<sub_number0;r++){
        fprintf(fp2, "%lf\t%lf\t%lf\n", Band0[r]);
    }
    
    fprintf(fp2, "%d\n", sub_number1);
    
    for(r=0;r<sub_number1;r++){
        fprintf(fp2, "%lf\t%lf\t%lf\n", Band1[r]);
    }
    */
    
    fclose(fp);
    //fclose(fp2);
    return;
}

//�ȉ����̂Ƃ���g��Ȃ��֐�

/*
void read_data(){
	    FILE *fp0;
    //int r_number0, r_number1;
    //int i;
    
    fp0 = fopen("savedata90.txt","r");
    
    fscanf(fp0, "%d", &r_number0);
    
//    double Band0[r_number0][3];
    
    for (i=0;i<r_number0;i++){
        fscanf(fp0,"%lf", &Band0[i][0]);
        fscanf(fp0,"%lf", &Band0[i][1]);
        fscanf(fp0,"%lf", &Band0[i][2]);
    }
    
    fscanf(fp0, "%d" ,&r_number1);
    
    r_number2 = r_number1;
    
//    double Band1[r_number1][3];
    
    for (i=0;i<r_number1;i++){
        fscanf(fp0,"%lf", &Band1[i][0]);
        fscanf(fp0,"%lf", &Band1[i][1]);
        fscanf(fp0,"%lf", &Band1[i][2]);
    }
    
    fclose(fp0);
    
    
    for(i=0;i<r_number0;i++){
        printf("Cr0=%d : %lf %lf %lf\n",i, Band0[i][0], Band0[i][1], Band0[i][2]);
    }
    for(r=1;r<r_number0-1;r++){
        printf("Cr0=%d : %f %f %f\n", r_number0-1+r, -Band0[r_number0-1-r][0], -Band0[r_number0-1-r][1], Band0[r_number0-1-r][2]);
    }

    printf("\n");
    
    for(i=0;i<r_number1;i++){
        printf("Cr1=%d : %lf %lf %lf\n", i,Band1[i][0], Band1[i][1], Band1[i][2]);
    }
    for(r=1;r<r_number1-1;r++){
        printf("Cr1=%d : %f %f %f\n", r_number1-1+r, -Band1[r_number1-1-r][0], -Band1[r_number1-1-r][1], Band1[r_number1-1-r][2]);
    }

	return;
}

int deltalength_matching(double band_delta_min, int r_number2,  double Band1[r_number1][3], double Band2[r_number2][3]){
    double b1tob2 = 0;
    double hosei[3];
    int r;
    for(r=0;r<r_number2;r++){
        printf("!!!!!!!Best point not found !!!!!!\n");
        b1tob2 = sqrt(pow((Band2[r][0]-Band1[r][0]), 2) + pow((Band2[r][1]-Band1[r][1]), 2) + pow((Band2[r][2]-Band1[r][2]), 2));
        hosei[0] = (Band2[r][0]-Band1[r][0])*(band_delta_min/b1tob2);
        hosei[1] = (Band2[r][1]-Band1[r][1])*(band_delta_min/b1tob2);
        hosei[2] = (Band2[r][2]-Band1[r][2])*(band_delta_min/b1tob2);
    
        Band2[r][0] = (Band1[r][0] + hosei[0]);
        Band2[r][1] = (Band1[r][1] + hosei[1]);
        Band2[r][2] = (Band1[r][2] + hosei[2]);
    }
    
    return 1;
}

*/