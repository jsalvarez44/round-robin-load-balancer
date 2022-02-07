#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <windows.h>
#include <string.h>

void saludo();
void presentacion();
int  interfazEntrada();
bool verificacionUsuario();
void salir();
void ingresarCaso();
void modificarCaso();
void buscarCaso();
void eliminarCaso();
void salida();

main(){
	saludo();
	presentacion();
	int eleccion;
	
	for(;;){
		system("cls");
		if(verificacionUsuario()){
			for(;;){
				system("cls");
				eleccion = interfazEntrada();
			
				switch(eleccion){
					case 1: ingresarCaso();break;
					case 2: modificarCaso();break;
					case 3: buscarCaso();break;
					case 4: eliminarCaso();break;
					case 5: salir();break;
					default: printf("\nIngrese una opcion correcta\n");getch();break;
				}
			}
		}
		else{
			salida();
		}
	}
	
}

void saludo(){
	printf("\n\t\t___________________________________________________________\n\n");
	printf ("\t\t\t\tUNIVERSIDAD DE LAS FUERZAS ARMADAS\n");
	printf ("\t\t**********\t *********\t********* \t***********\n");
    printf ("\t\t**********\t**********\t**********\t***********\n");
    printf ("\t\t***       \t***       \t***    ***\t***\n");
    printf ("\t\t***       \t***       \t***    ***\t***\n");
    printf ("\t\t********  \t********* \t**********\t*********  \n");
    printf ("\t\t********  \t *********\t********* \t*********  \n");
    printf ("\t\t***       \t       ***\t****      \t***\n");
    printf ("\t\t***       \t       ***\t****      \t***\n");
    printf ("\t\t**********\t**********\t****      \t***********\n");
    printf ("\t\t**********\t********* \t****      \t***********\n"); 
    printf("\n\t\t___________________________________________________________");
    getch();
    system("cls");
}

void presentacion(){
	printf("\n\t\t______________________________________________\n\n");
    printf("\t\t**********    \t**********    \t***    *** \n");
    printf("\t\t**********    \t**********    \t***    *** \n");
    printf("\t\t   ****       \t   ****       \t***    *** \n");
    printf("\t\t   ****       \t   ****       \t***    *** \n");
    printf("\t\t   ****       \t   ****       \t***    *** \n");
    printf("\t\t   ****       \t   ****       \t***    *** \n");
    printf("\t\t   ****    ***\t   ****    ***\t********** ***\n");
    printf("\t\t   ****    ***\t   ****    ***\t********** ***\n");
    printf("\t\t\t\t ABOGADOS");
    printf("\n\t\t______________________________________________");
    getch();
}

int interfazEntrada(){
	int opc;
	
	printf("\n\tMENU DE INICIO");
	printf("\n\tIngrese el numero de su eleccion");
	printf("\n\t__________________________________");
	printf("\n\n\t1. Ingresar Caso");
	printf("\n\t2. Modificar Caso");
	printf("\n\t3. Buscar Caso");
	printf("\n\t4. Eliminar Caso");
	printf("\n\t5. Salir");
	printf("\n\t__________________________________\n\n");
	
	printf("\tOpcion: ");
	scanf("%d",&opc);
	return opc;
}

bool verificacionUsuario(){
	char usuario[100],contrasena[100];
	printf("\n\tVERIFICACION DE USUARIO\n");
	printf("\t_______________________________\n\n");
	printf("\tIngrese usuario: ");
	scanf("%s",usuario);
	printf("\tIngrese contraseña: ");
	scanf("%s",contrasena);
	printf("\t_______________________________");
	
	if(strcmp(usuario,"Admin")==0 && strcmp(contrasena,"Passw")==0){
		printf("\n\n\tContraseña correcta...");
		Sleep(3000);
		return true;
	}
	else{
		return false;
	}
}

void ingresarCaso(){
	printf("\tIngresando Caso...");
	getch();
}

void modificarCaso(){
	printf("\tModificando Caso...");
	getch();
}

void buscarCaso(){
	printf("\tBuscando Caso...");
	getch();
}

void eliminarCaso(){
	printf("\tEliminando Caso...");
	getch();
}

void salir (void){
	printf("\nFinalizando programa...\n");
	Sleep(3000);
	exit(1);
}

void salida(){
	int opc;
	printf("\n\n\tContraseña Incorrecta\n");
	printf("\tDesea salir del sistema? \n\n\t1.Si\n\t2.No\n\n");
	printf("\tOpcion: ");
	scanf("%d",&opc);
	
	if(opc == 1){
		salir();
	}
	else if(opc == 2){
		return;
	}
}


















