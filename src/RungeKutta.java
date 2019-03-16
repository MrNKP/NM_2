import java.io.*;

import static java.lang.Math.abs;
//import com.lab2.io.FileWrite;

/*
    Задание:
        Решение задачи Коши с заданной точностью с автоматическим
        выбором шага методом удвоения и деления шага пополам.

        Интегрирование обыкновенного дифференциального уравнения
            y' = f(x,y), x in [А,В]
        с начальным условием y(c)=yc,
        где точка cовпадает либо с началом,
        либо с концом отрезка интегрирования.

        Состав входного файла:
        1. A B x y
        2. hmin epsmax
 */
class RungeKutta {
    private FileWrite fileOut;  // вывод данных в файл
    private static final double mconst = 0.25; // 1/2^2
    private int error;          // индикатор ошибки: 0 - ок, 1 - ошибка ввода
    private boolean direction;  // false - left-to-right, true - right-to-left
    private double h;           // текущий шаг интегрирования
    private double
            A, B,
            x, y,   // c ~ x, yc ~ y
            hmin, hmax,   // минимальный шаг
            eps,    // оценка погрешности на шаге
            epsmax; // наибольшая допустимая погрешность
    private int pointsCount, wAccuracy; // кол-во точек, кол-во точек в которых не достигли точности

    // вывод данных в консоль
    private <T> void write(T ...items) {
        for (T i : items)
            System.out.print(i + "\t");
    }
    private <T> void writeln(T ...items) {
        for (T i : items)
            System.out.print(i + "\t");
        System.out.println();
    }

    private void dataFromFile(String inputfile, String outfile) {
        fileOut = new FileWrite(outfile);
        // считывание из файла. try-with-resources конструкция, закрывает поток сам
        try(BufferedReader fin = new BufferedReader(new FileReader(new File(inputfile)))){
            A = Double.parseDouble(fin.readLine());
            B = Double.parseDouble(fin.readLine());
            x = Double.parseDouble(fin.readLine());
            y = Double.parseDouble(fin.readLine());
            hmin = Double.parseDouble(fin.readLine());
            epsmax = Double.parseDouble(fin.readLine());
        }
        catch(IOException ex){
            //System.out.printf("File error!");
            System.out.println(ex.getMessage());
            ex.printStackTrace();
        }

    }

    public RungeKutta() {
        this(
                "d://IdeaProjects//NM_2//data.txt",
                "d://IdeaProjects//NM_2//output.txt"
        );
    }
    public RungeKutta(String inputfile, String outfile) {
        wAccuracy = error = 0;
        pointsCount = 1;
        dataFromFile(inputfile, outfile);
        hmax = h = abs(B-A)/10.0;
        direction = (x == B); // определение направления

        if ((direction ? B : A) != x || A >= B || epsmax < 0 || hmin <= 0)
            error = 1;
    }

    // примеры функций
    private static class functs{
        static double f1(double X, double Y) {
            return 12*X; // 6*x^2
        }
        static double f2(double X, double Y) {
            return 1; // X
        }
        static double f3(double X, double Y) {
            return 6; // X
        }
    }
    private double f(double X, double Y) {
        return functs.f1(X, Y);    // ф-ия изменяется при необходимости
    }

    // в зависимости от направления выдает -1(если справа налево) или 1(слева направо)
    private double direct() {
        return (direction) ? -1.0 : 1.0;
    }

    // Нахождение yi+1. Вычисление по формуле (110), метод 3-го порядка
    private double intMethod113(double X, double Y, double H) {
        // TODO добавить вычисление при right-to-left direction
        double K1 = H*f(X, Y);
        double K2 = H*f(X + (1.0/2.0)*H, Y + (1.0/2.0)*K1);
        return Y + K2*direct();
    }

    // Оценка погрешности по правилу Рунге
    private double rungeRule(double yh1, double yh2) {
        return abs((yh1 - yh2)/(mconst - 1));
    }

    // Вычисление Yi+1 + оценка погрешности
    private double runge() {
        // 1) вычисление Yi+1
        double nextY = intMethod113(x, y, h); // вычисление yi+1

        // 2) Оценка погрешности по правилу Рунге
        // вычисление первого h/2 шага внутри второго шага
        eps = rungeRule(nextY, intMethod113(x + h/2, intMethod113(x, y, h/2), h/2));
        //if (eps < 1e-14) eps = 0;
        return nextY;
    }

    private void step() {
        y = runge();        // вычисление Yi+1, оценка погрешности по правилу Рунге
        if (eps > epsmax)
            wAccuracy++;
        x += h * direct();  // делаем шаг с учётом направления
        fileOut.write(3, x, y, eps, h); // выводим данные
    }
    private void laststep() {
        h = direction ? x - A : B - x;
        step();
    }



    public void start() {
        if (error == 1) {
            System.out.println("Ошибка ввода!");
            return;
        }

        boolean lastdiv = false; // было ли последним действием деление шага

        //h = 0.81;    // искуственный шаг
        fileOut.cleanFile();
        fileOut.write("X\t\t\tY\t\t\tACC\t\t\tH");
        fileOut.write(3, x, y, eps); // пишем первое граничное значение
        // Основной цикл
        while (true) {

            // автовыбор шага
            runge();
            if (eps <= epsmax/8)     // если точность слишком высокая
                while (eps <= epsmax/8)
                    // наибольший шаг взят условно, что бы сильно не "уплыло"
                    // однако, возможно, это ограничение вовсе не нужно
                    if (h * 2 <= hmax) {
                        if (lastdiv) { // не допускаем умножения после деления
                            lastdiv = false;
                            break;
                        }
                        h *= 2;
                        runge();
                        // если при умножении вышли за макс погрешность, то откатываемся
                        if (eps > epsmax) {
                            h /= 2;
                            break;
                        }
                    }
                    else
                        break;
            else
            if (eps > epsmax)   // если точность слишком низкая
                while (eps > epsmax)
                    if (h / 2 >= hmin) {
                        lastdiv = true;
                        h /= 2;
                        runge();
                    }
                    else { // если не удалось достичь точности
                        h = hmin; // устанавливаем мин шаг
                        //wAccuracy++; // прибавляем кол-во точек в которых не достигли точности
                        break;
                    }

            // Определение конца отрезка.
            if ((direction ? x-h-A : B-x-h) < hmin || x + h == x) // Если после текущего шага B-x будет < hmin
                break;

            /*
            y = runge();
            x += h * direct();  // делаем шаг с учётом направления
            fileOut.write(3, x, y, eps, h); // выводим данные
            */
            step();

            pointsCount++;
        }

        // возможные проблемы с определением погрешности на границе
        // могут быть связаны с сильно дробным шагом

        // обрабатываем состояние, когда находимся у конца отрезка
        // выбор направления - direction ? [справа налево] : [слева направо]
        if ((direction ? x-A : B-x) >= 2*hmin) { // если больше чем за 2 мин шага от края
            h = direction ? x - hmin - A : B - hmin - x;
            step();
            laststep();
            pointsCount++;
        }
        else
        if ((direction ? x-A : B-x ) <= 1.5*hmin) // если меньше чем за 1.5 мин шага
            laststep();
        else { // если 1.5*hmin < остаток(B-x) < 2*hmin
            h = direction ? (x-A)/2.0 : (B-x)/2.0;
            step();
            laststep();
            pointsCount++;
        }
        pointsCount++;

        fileOut.write("\n");
        fileOut.write("Points count:   \t" + pointsCount);
        fileOut.write("Accuracy bad in: \t" + wAccuracy);
    }
}

