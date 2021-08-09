function action(input, output, filename) {
        open(input + filename);
        makeRectangle(460, 136, 670, 700);
        run("Crop");
        saveAs("tiff", output + filename);
        close();
}

path1 = "C:/Users/15183/Documents/Spring_data_2021/8-7-21_1388_Victoria/control_200 (1)";
path2 = "C:/Users/15183/Documents/Spring_data_2021/images/Images";
input = File.getDirectory(path1);
output = File.getDirectory(path2);

setBatchMode(true); 
list = getFileList(input);
for (i = 0; i < list.length; i++){
        action(input, output, list[i]);
}
setBatchMode(false);



## if you want to know what files are in your directory folder (make sure you're in the right place)
path1 = "C:/Users/15183/Documents/Spring_data_2021/8-7-21_1388_Victoria/control_200 (1)";
input = File.getDirectory(path1);

setBatchMode(true); 
list = getFileList(input);
for (i = 0; i < list.length; i++){
        print(list[i]);
}
setBatchMode(false);

