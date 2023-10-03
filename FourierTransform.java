import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import java.util.List;
import java.util.stream.*;

public class FourierTransform {

    // Часовий крок (інтервал) для побудови графіку.
    private static final int T = 5;

    public static void main(String[] args) throws IOException {
        System.out.println("Читаємо дані з файлу...");

        double[] Y;
        double[] X;
        // Зчитуємо дані з файлу f13.txt, розділяємо їх за пробілами та перетворюємо у масив double
        try (Stream<String> stream = Files.lines(Paths.get("f13.txt"))) {
            Y = stream.flatMap(line -> Arrays.stream(line.split("\\s+")))
                    .mapToDouble(Double::parseDouble)
                    .toArray();
            X = IntStream.range(0, Y.length)
                    .mapToDouble(i -> T * ((double) i / Y.length))
                    .toArray();
        }

        System.out.println("Дані завантажено. Кількість точок: " + Y.length);
        System.out.println("Будуемо графік \"Початкові дані\"...");
        plotGraph("Початкові дані", X, Y);

        System.out.println("Знаходимо найближчий ступінь двійки...");
        // Визначаємо оптимальний розмір масиву для перетворення Фур'є (ступінь 2)
        int paddedLength = 1;
        while (paddedLength < Y.length) {
            paddedLength *= 2;
        }
        System.out.println("Найближчий ступінь двійки: " + paddedLength);


        System.out.println("Додаємо нулі до масиву даних...");
        double[] paddedY = new double[paddedLength];
        System.arraycopy(Y, 0, paddedY, 0, Y.length);

        System.out.println("Виконуємо перетворення Фур'є...");
        // Виконуємо швидке перетворення Фур'є та отримуємо модуль результатів
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] complexResults = fft.transform(paddedY, TransformType.FORWARD);
        double[] Fur = Arrays.stream(complexResults).mapToDouble(Complex::abs).toArray();
        System.out.println("Будуемо графік \"Перетворення Фур'є\"...");
        plotGraph("Перетворення Фур'є", X, Fur);

        System.out.println("Беремо половину результатів перетворення Фур'є...");
        double[] Furh = Arrays.copyOfRange(Fur, 0, Y.length / 2);
        double[] N1 = IntStream.range(0, Furh.length).mapToDouble(i -> (double) i / T).toArray();
        double[] N2 = Arrays.copyOfRange(N1, 0, Y.length / 2);
        System.out.println("Будуемо графік \"Половина перетворення Фур'є\"...");
        plotGraph("Половина перетворення Фур'є", N2, Furh);

        System.out.println("Пошук локальних максимумів...");
        // Визначаємо локальні максимуми в результаті Фур'є
        List<Integer> extr = new ArrayList<>();
        for (int i = 9; i < Furh.length; i++) {
            boolean isMax = true;
            for (int j = 1; j <= 9; j++) {
                if (!(Furh[i] > Furh[i - j])) {
                    isMax = false;
                    break;
                }
            }
            if (isMax) {
                extr.add(i);
            }
        }
        extr = extr.stream().filter(e -> e >= T).collect(Collectors.toList());
        int max1 = Math.round((float)extr.get(0) / T);
        System.out.println("Знаходимо компоненти найбільшого впливу. Існує в точці:");
        System.out.println("𝜗1 = " + max1);

        System.out.println("Будуємо та розв’язуємо систему рівнянь, щоб знайти коефіцієнти:");
        //Будуємо та розв'язуємо систему рівнянь для апроксимації даних
        double[][] coefficientsArray = new double[5][5];
        double[] constantsArray = new double[5];

        double[][][] computedValues = {
                // Отримуємо коефіцієнти для рівнянь
                compute(X, Y, max1, 1),
                compute(X, Y, max1, 2),
                compute(X, Y, max1, 3),
                compute(X, Y, max1, 4),
                compute(X, Y, max1, 5)
        };

        for (int i = 0; i < 5; i++) {
            coefficientsArray[i] = computedValues[i][0];
            constantsArray[i] = computedValues[i][1][0];
        }

        String[] equations = {"x^3", "x^2", "x", "sin(2π" + max1 + "x)", "const"};
        for (int i = 0; i < 5; i++) {
            StringBuilder equation = new StringBuilder();
            for (int j = 0; j < 5; j++) {
                equation.append(coefficientsArray[j][i]).append(" * ").append(equations[j]).append(" + ");
            }
            equation.delete(equation.length() - 3, equation.length());
            equation.append("= ").append(constantsArray[i]);
            System.out.println(equation);
        }

        RealMatrix coefficients = new Array2DRowRealMatrix(5, 5);
        double[] constants = new double[5];

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                coefficients.setEntry(i, j, computedValues[j][0][i]);
            }
            constants[i] = computedValues[i][1][0];
        }

        // Використовуємо метод LUDecomposition для розв'язання системи лінійних рівнянь
        DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
        RealVector solution = solver.solve(new ArrayRealVector(constants, false));

        System.out.println("Розв'язуємо рівняння");
        for (int i = 0; i < 5; i++) {
            System.out.println("Розв'язок " + (i + 1) + ": " + solution.getEntry(i));
        }

        // Обчислюємо апроксимовані значення y за допомогою знайдених коефіцієнтів
        double[] yApproximated = new double[X.length];
        for (int i = 0; i < X.length; i++) {
            yApproximated[i] = solution.getEntry(0) * Math.pow(X[i], 3)
                    + solution.getEntry(1) * Math.pow(X[i], 2)
                    + solution.getEntry(2) * X[i]
                    + solution.getEntry(3) * Math.sin(2 * Math.PI * max1 * X[i])
                    + solution.getEntry(4);
        }
        // Виводимо отриману апроксимуючу функцію
        System.out.println("Отримана апроксимуюча функція:");
        String function = solution.getEntry(0) + " * x^3 + "
                + solution.getEntry(1) + " * x^2 + "
                + solution.getEntry(2) + " * x + "
                + solution.getEntry(3) + " * sin(2π" + max1 + "x) + "
                + solution.getEntry(4);
        System.out.println("y(x) = " + function);

        double error = 0.0;
        for (int i = 0; i < Y.length; i++) {
            error += Math.pow(yApproximated[i] - Y[i], 2);
        }
        System.out.println("Помилка: " + error);

        System.out.println("Будуемо графік \"Наближена функція\"...");
        plotGraph("Наближена функція", X, yApproximated);
    }

    /**
     * Обчислює коефіцієнти для рівнянь системи, які використовуються для апроксимації даних.
     *
     * @param X     масив значень X з вхідних даних.
     * @param Y     масив значень Y з вхідних даних.
     * @param max1  основна частота для перетворення Фур'є.
     * @param index індекс, що визначає тип рівняння (від 1 до 5).
     * @return      двомірний масив, де перший підмасив включає коефіцієнти a1, a2, a3, a4, a5,
     *              а другий підмасив містить одне значення, що відповідає правій частині рівняння.
     */
    public static double[][] compute(double[] X, double[] Y, int max1, int index) {
        double a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, minus = 0;
        for (int i = 0; i < X.length; i++) {
            double x = X[i];
            double y = Y[i];

            switch (index) {
                case 1 -> {
                    a1 += Math.pow(x, 6);
                    a2 += Math.pow(x, 5);
                    a3 += Math.pow(x, 4);
                    a4 += Math.pow(x, 3) * Math.sin(2 * Math.PI * max1 * x);
                    a5 += Math.pow(x, 3);
                    minus += Math.pow(x, 3) * y;
                }
                case 2 -> {
                    a1 += Math.pow(x, 5);
                    a2 += Math.pow(x, 4);
                    a3 += Math.pow(x, 3);
                    a4 += Math.pow(x, 2) * Math.sin(2 * Math.PI * max1 * x);
                    a5 += Math.pow(x, 2);
                    minus += Math.pow(x, 2) * y;
                }
                case 3 -> {
                    a1 += Math.pow(x, 4);
                    a2 += Math.pow(x, 3);
                    a3 += Math.pow(x, 2);
                    a4 += x * Math.sin(2 * Math.PI * max1 * x);
                    a5 += x;
                    minus += x * y;
                }
                case 4 -> {
                    a1 += Math.pow(x, 3) * Math.sin(2 * Math.PI * max1 * x);
                    a2 += Math.pow(x, 2) * Math.sin(2 * Math.PI * max1 * x);
                    a3 += x * Math.sin(2 * Math.PI * max1 * x);
                    a4 += Math.sin(2 * Math.PI * max1 * x) * Math.sin(2 * Math.PI * max1 * x);
                    a5 += Math.sin(2 * Math.PI * max1 * x);
                    minus += Math.sin(2 * Math.PI * max1 * x) * y;
                }
                case 5 -> {
                    a1 += Math.pow(x, 3);
                    a2 += Math.pow(x, 2);
                    a3 += x;
                    a4 += Math.sin(2 * Math.PI * max1 * x);
                    a5 += 1;
                    minus += y;
                }
            }
        }
        return new double[][]{{a1, a2, a3, a4, a5}, {minus}};
    }

    // Допоміжний метод для побудови графіків з використанням JFreeChart
    public static void plotGraph(String title, double[] xData, double[] yData) {
        XYSeries series = new XYSeries("Дані");
        for (int i = 0; i < xData.length; i++) {
            series.add(xData[i], yData[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                title,
                "X",
                "Y",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();

        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        Dimension frameSize = frame.getSize();
        int x = (screenSize.width - frameSize.width) / 2;
        int y = (screenSize.height - frameSize.height) / 2;
        frame.setLocation(x, y);

        frame.setVisible(true);
    }

}


