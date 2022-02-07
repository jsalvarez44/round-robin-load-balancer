#include<stdio.h>
#include<stdlib.h>
#include <conio.h>
#include <string.h>

char opc,opc1;
void presentacion();
void ingreso();
void contra ();
void menuprincipal(char opc,char opc1);
void menu2(char opc1);

int main (){
    char opc,opc1;
    presentacion();
    ingreso();
    return 0;
}
void menu2(char opc1){
    do{
        printf("\t***MENU***\na.-Personas naturales\nb.-Personas juridicas\nc.-Salir\nSu opcion es:");
        scanf("%s",&opc1);
        switch (opc1){
            case 'a':
                printf("demandados\n");
            break;
            case 'b':
                printf("tramites\n");
            break;
        }
    }while (opc1!='c');
    system("cls");
}
void menuprincipal(char opc,char opc1){
    do{
        printf("\t***MENU PRINCIPAL***\na.-Ingresar caso\nb.-Modificar caso\nc.-Buscar caso\nd.-Borrar caso\ne.-salir\nIngrese su opcion:");
        scanf("%s",&opc);
        switch (opc){
        case 'a':
            menu2(opc1);
            break;
        case 'b':
            menu2(opc1);
            break;
        case 'c':
            menu2(opc1);
            break;
        case 'd':
            menu2(opc1);
            break;
        }
    }while (opc!='e');
    system("cls");
}
void contra ( ){
char pass[20];
    printf("\nIngrese su Clave: ");
    int conta=0;
    int i=0;

     //Se piden caracteres hasta que se introduca enter
    while(pass[i]!=13)
    {

         //Capturamos carácter
        pass[i]=getch();

        //Si es un carácter válido y no se ha sobrepasado el límite de 20 caracteres se imprime un asterisco
        if(pass[i]>32 && i<20)
                {
                  putchar('*');
                  i++;
                }
        //Si se pulsa la tecla RETROCESO, se retrocede un carácter, se imprime un espacio para eliminar el asterisco y se vuelve a retroceder para que el siguiente asterisco se coloque a continuación del anterior.
        else if(pass[i]==8 && i>0)
                {
                  putchar(8);
                  putchar(' ');
                  putchar(8);
                  i--;
                }
    }

    pass[i]='\0';

    if(!strcmp(pass, "loquesea")){
        printf("\n\n\n\n\n\t\tACCESO PERMITIDO\n");
        menuprincipal( opc, opc1);
    }else{
        printf("\n\n\n\n\n\t\t**ACCESO DENEGADO**");
    }


getch();
}
void presentacion(){
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
     printf ("\t\t\t\tUNIVERSIDAD DE LAS FUERZAS ARMADAS\n");
    system ("pause");
    system ("cls");
    printf("\t\t**********    \t**********    \t***    *** \n");
    printf("\t\t**********    \t**********    \t***    *** \n");
    printf("\t\t   ****       \t   ****       \t***    *** \n");
    printf("\t\t   ****       \t   ****       \t***    *** \n");
    printf("\t\t   ****       \t   ****       \t***    *** \n");
    printf("\t\t   ****       \t   ****       \t***    *** \n");
    printf("\t\t   ****    ***\t   ****    ***\t********** ***\n");
    printf("\t\t   ****    ***\t   ****    ***\t********** ***\n");
    printf("\t\t\t\t ABOGADOS \n\n\n");
    system ("pause");
    system ("cls");
}
void ingreso(){
    char auxttu;
    printf("DESEA CONTINUAR s/n?:");
    scanf ("%s",&auxttu);
    if (auxttu=='s'){
        contra();
    }
    else if (auxttu=='n'){
        printf("ADIOS =(");
    }else{
        printf ("OPCION INCORRECTA");
    }
}
