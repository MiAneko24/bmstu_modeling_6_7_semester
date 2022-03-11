import kotlin.math.abs
import kotlin.math.ceil
import kotlin.math.pow
import jetbrains.letsPlot.export.ggsave
import jetbrains.letsPlot.geom.geomDensity
import jetbrains.letsPlot.geom.geom_density
import jetbrains.letsPlot.ggsize
import jetbrains.letsPlot.letsPlot
import jetbrains.letsPlot.lets_plot
import javax.swing.*
import org.math.plot.Plot2DPanel
import java.awt.Color

val output_step = 1

fun picarFirstApprox(x : Double, y0 : Double) : Double = y0 + x.pow(3) / 3 + y0 * y0

fun picarSecondApprox(x : Double, y0 : Double) : Double = y0 + x.pow(3) / 3 + x.pow(7) / 63

fun picarThirdApprox(x : Double, y0 : Double) : Double = y0 + x.pow(3) / 3 + x.pow(7) / 63 + 2 * x.pow(11) /
        (21 * 9 * 11) + x.pow(15) / (63 * 63 * 15)

fun picarFourthApprox(x : Double, y0 : Double) : Double = y0 + x.pow(3) / 3 + x.pow(7) / 63 + 2 * x.pow(11) / (3 * 63 * 11) +
        117 * x.pow(15) / (63 * 63 * 3 * 11 * 15) + (41 * 2 * x.pow(19)) / (63 * 63 * 3 * 11 * 15 * 19) +
        (2 * 11 * 11 * 3 + 4 * 63 * 5) * x.pow(23) / (63.0.pow(3) * 9 * 5 * 11 * 11 * 23) + 4 * x.pow(27) /
        (3 * 63.0.pow(3) * 11 * 15 * 27) + x.pow(31) / (63.0.pow(4) * 15 * 15 * 31)

fun picardSolver(x0 : Double, y0 : Double, h : Double, n : Int) : MutableList<MutableList<Double>> {
    val resList : MutableList<MutableList<Double>> = mutableListOf()
    var x = x0
    var cnt = 0
    for (i in 0..n) {
        if (cnt % output_step == 0)
            resList.add(mutableListOf(x, picarFirstApprox(x, y0), picarSecondApprox(x, y0), picarThirdApprox(x, y0), picarFourthApprox(x, y0)))
        x += h
        cnt++
    }
    return resList
}

fun f(x : Double, y : Double) : Double = x * x + y * y


fun eulerSolver (x0: Double, y0 : Double, h : Double, n : Int) : MutableList<Double> {
    val resList : MutableList<Double> = mutableListOf()

    var x = x0
    var y = y0
    var cnt = 0
    for (i in 0..n) {
        if (cnt % output_step == 0)
            resList.add(y)
        y += h * f(x, y)
        x += h
        cnt++
    }
    return resList
}

fun rungeKuttSolver(x0 : Double, y0 : Double, alpha : Double, h : Double, n : Int) : MutableList<Double> {
    val resList : MutableList<Double> = mutableListOf()
    var x = x0
    var y = y0
    var cnt = 0
    for (i in 0..n) {
        if (cnt % output_step == 0)
            resList.add(y)
        y += h * ((1 - alpha) * f(x, y) + alpha * f(x + h / (2 * alpha), y + h / (2 * alpha) * f(x, y)))
        x += h
        cnt++
    }
    return resList
}

fun printResult(picar : MutableList<MutableList<Double>>, euler : MutableList<Double>, runge : MutableList<Double>) {
    println("╔═════════════════╦═════════════════╦═════════════════╦═════════════════╦═════════════════╦═════════════════╦═══════════════════╗")
    println("║         X       ║   Метод Пикара  ║   Метод Пикара  ║   Метод Пикара  ║   Метод Пикара  ║   Метод Эйлера  ║ Метод Рунге-Кутта ║")
    println("║                 ║(1ое приближение)║(2ое приближение)║(3-е приближение)║(4ое приближение)║                 ║                   ║")
    println("╠═════════════════╬═════════════════╬═════════════════╬═════════════════╬═════════════════╬═════════════════╬═══════════════════╣")
    for (i in 0 until picar.size) {
        println("║    %9.4f    ║    %9.4f    ║    %9.4f    ║    %9.4f    ║    %9.4f    ║    %9.4f    ║     %9.4f     ║".format(
            picar[i][0], picar[i][1], picar[i][2], picar[i][3], picar[i][4], euler[i], runge[i]))
        if (i != picar.size - 1) {
            println("╠═════════════════╬═════════════════╬═════════════════╬═════════════════╬═════════════════╬═════════════════╬═══════════════════╣")
        }
        else
            println("╚═════════════════╩═════════════════╩═════════════════╩═════════════════╩═════════════════╩═════════════════╩═══════════════════╝")
    }
}

fun main() {
    val xStart = 0.0
    val xEnd = 2.0
    val y0 = 0.0
    val h = 10e-3

    val n : Int = ceil(abs(xEnd - xStart) / h).toInt()

    val picar = picardSolver(xStart, y0, h, n)
    val euler = eulerSolver(xStart, y0, h, n)
    val runge = rungeKuttSolver(xStart, y0, 0.5, h, n)
    printResult(picar, euler, runge)
    val xgraph = mutableListOf<Double>()
    val npicar = picardSolver(xStart, y0, -h, n)
    val neuler : MutableList<Double> = eulerSolver(xStart, y0, -h, n).asReversed()
    val nrunge = rungeKuttSolver(xStart, y0, 0.5, -h, n).asReversed()
//    printResult(npicar, neuler, nrunge)
    val gpic = mutableListOf(mutableListOf(), mutableListOf(), mutableListOf(), mutableListOf<Double>())
    for (i in npicar.reversed())
    {
        xgraph.add(i[0])
        gpic[0].add(i[1])
        gpic[1].add(i[2])
        gpic[2].add(i[3])
        gpic[3].add(i[4])
    }

    for (i in picar)
    {
        xgraph.add(i[0])
        gpic[0].add(i[1])
        gpic[1].add(i[2])
        gpic[2].add(i[3])
        gpic[3].add(i[4])
    }

    neuler.addAll(euler)
    nrunge.addAll(runge)
    print("xl = ${xgraph.size}, gp[i]l = ${gpic[0].size}, el = ${neuler.size}, rl = ${nrunge.size}")
    val plot = Plot2DPanel()
    plot.addLinePlot("Picar 1st approx", Color.CYAN, xgraph.toDoubleArray(), gpic[0].toDoubleArray())
    plot.addLinePlot("Picar 2nd approx", Color.GREEN, xgraph.toDoubleArray(), gpic[1].toDoubleArray())
    plot.addLinePlot("picar 3rd approx", Color.ORANGE, xgraph.toDoubleArray(), gpic[2].toDoubleArray())
    plot.addLinePlot("Picar 4th approx", Color.MAGENTA, xgraph.toDoubleArray(), gpic[3].toDoubleArray())
    plot.addLinePlot("Euler", Color.RED, xgraph.toDoubleArray(), neuler.toDoubleArray())
    plot.addLinePlot("Runge-Kutt", Color.YELLOW, xgraph.toDoubleArray(), nrunge.toDoubleArray())
    val frame = JFrame()
    plot.addLegend("SOUTH")

    frame.setSize(600, 600)
    frame.contentPane = plot
    frame.isVisible = true

    frame.defaultCloseOperation = 2

}