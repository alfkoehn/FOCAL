#include "save_data.h"

void create_folder(saveData *saveDCfg){

    simulation_folder( projectPath );
    data_folder( projectPath, foldername );
    copyJSON( projectPath, foldername );

    exit(-1);
    
}

void simulation_folder(const char *path){
    
    struct stat st = {0};

    /*Checks if directory exists*/
    if( stat(path, &st) == -1){
        //Directory does not exists. Create it.
        if(mkdir(path, 0700) == 0){
            printf("Main project folder created successfully.\n");
        }else{
            printf("Error creating directory: %s\n", path);
            return;
        }
    }
}

void data_folder(const char *path, const char *folder_name){
    
    char fullPath[1024];

    // Create the full directory path and check for buffer overflow
    if (snprintf(fullPath, sizeof(fullPath), "%s/%s", path, folder_name) >= sizeof(fullPath)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
        return;
    }

    struct stat st = {0};

    /*Checks if directory exists.*/
    if( stat(fullPath, &st) == -1){
        //Directory does not exists. Create it.
        if( mkdir(fullPath, 0700) == 0){
            printf("%s folder created successfully. \n", folder_name);
        }else{
            printf("Error creating directory: %s\n", folder_name);
            return;
        }
    }else{
        printf("%s already exists.\n", folder_name);
    }
}

void copyJSON(const char *path, const char *folder_name){
    
    char destination[1024];

    //Read the source file
    FILE *srcFile = fopen("input_FOCAL.json", "rb");
    if (srcFile == NULL) {
        perror("Error opening source file");
        return;
    }

    // Create the full directory path and check for buffer overflow
    if (snprintf(destination, sizeof(destination), "%s/%s/input_FOCAL.json", path, folder_name) >= sizeof(destination)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
        return;
    }

    //open destination file
    FILE *destFile = fopen(destination, "wb");
    if (destFile == NULL) {
        perror("Error opening destination file");
        fclose(srcFile);
        return;
    }

    char buffer[1024];
    size_t bytesRead;
    while ((bytesRead = fread(buffer, 1, sizeof(buffer), srcFile)) > 0) {
        fwrite(buffer, 1, bytesRead, destFile);
    }

    fclose(srcFile);
    fclose(destFile);

    printf("JSON file saved.\n");
}

