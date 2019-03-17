import java.io.*;

public class FileWrite {
    // pathname of output files
    private String filepath;
    PrintWriter pw = new PrintWriter(System.out, true);

    public FileWrite(String path) {
        filepath = path;
    }

    public void write(int precision, double... tmp) {
        File fout = new File(filepath);
        try(PrintWriter file = new PrintWriter( new FileWriter(filepath, true))){
            String template = "%1$." + Integer.toString(precision) + "f";

            if(tmp.length == 0 || tmp.length > 4)
                file.print("Неправильное кол-во аргументов: " + tmp.length);
            else {
                // ATTENTION костыли
                file.println();

                file.printf(template,tmp[0]);
                file.print("\t");
                if (tmp[0] < 1000)
                    file.print("\t");

                file.printf(template,tmp[1]);
                file.print("\t");
                if (tmp[1] < 1000)
                    file.print("\t");

                file.print(tmp[2]);

                if(tmp.length == 4){
                    file.print("\t");
                    if (tmp[3] < 1000)
                        file.print("\t");
                    file.printf(template,tmp[3]);
                }
            }
            //file.println();
        }catch (IOException exp){
            pw.println("Error!!!");
            pw.println("File writing handling problem!");
            exp.printStackTrace();
        }
    }

    public void write(double tmp){
        File fout = new File(filepath);
        try(BufferedWriter file = new BufferedWriter(new FileWriter(fout, true))){
            file.append(Double.toString(tmp));
            file.flush();
            file.newLine();
        }catch (IOException exp){
            pw.println("Error!!!");
            pw.println("File writing handling problem!");
            exp.printStackTrace();
        }
    }

    public void write(String str){
        File fout = new File(filepath);
        try(BufferedWriter file = new BufferedWriter(new FileWriter(fout, true))){
            file.append(str);
            file.flush();
            file.newLine();
        }catch (IOException exp){
            pw.println("Error!!!");
            pw.println("File writing handling problem!");
            exp.printStackTrace();
        }
    }

    public void cleanFile(){
        try(FileWriter fileWriter = new FileWriter(filepath)){
            fileWriter.write("");
        }catch (IOException exp){
            pw.println(exp.getStackTrace());
        }
    }
}