import org.math.plot.Plot2DPanel
import java.awt.*
import javax.swing.JFrame
import javax.swing.JPanel
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.pow
import kotlin.math.sign

// в зависимостях 
// implementation("com.github.yannrichet:JMathPlot:1.0.1")

class RungeKuttSolver(val R : Double, val m : Double, val k0 : Double, val T0 : Double, val Tw : Double, val p : Int) {
    val c = 3e10

    fun T(z : Double) : Double = (Tw - T0) * z.pow(p) + T0

    fun u_p(z : Double) : Double = 3.084e-4 / (exp(4.799e4 / T(z)) - 1)

    fun k(z : Double) : Double = k0 * (T(z) / 300).pow(2)

    fun u(z : Double, F: Double) = -F*3*R*k(z) / c

    fun F(z : Double, F : Double, u:Double) : Double =
        if (abs(z) < 1e-7)
            c * R / 2 * k(z) * (u_p(z) - u)
        else
            R * c * k(z) * (u_p(z) - u) - F / z

    fun runge_4th(z_start : Double, z_end : Double, n : Double, f_start : Double, u_start : Double) : MutableList<MutableList<Double>> {
        val res = mutableListOf<MutableList<Double>>(mutableListOf(z_start), mutableListOf(f_start), mutableListOf(u_start))
        var z_cur = z_start
        var f_cur = f_start
        var u_cur = u_start
        val h = 1 / (n - 1)
        while (abs(z_cur - z_end) > 1e-7) { //k - for f, q - for u
            val k1 = h * F(z_cur, f_cur, u_cur)
            val q1 = h * u(z_cur, f_cur)

            val k2 = h * F(z_cur + h/2, f_cur + k1/2, u_cur + q1/2)
            val q2 = h * u(z_cur + h/2, f_cur + k1/2)

            val k3 = h * F(z_cur + h/2, f_cur + k2/2, u_cur + q2/2)
            val q3 = h * u(z_cur + h/2, f_cur + k2/2)

            val k4 = h * F(z_cur + h, f_cur + k3, u_cur + q3)
            val q4 = h * u(z_cur + h, f_cur + k3)
            f_cur += (k1 + 2 * k2 + 2 * k3 + k4) / 6
            u_cur += (q1 + 2 * q1 + 2 * q3 + q4) / 6
            z_cur += h
            res[0].add(z_cur)
            res[1].add(f_cur)
            res[2].add(u_cur)
        }
        return res
    }

    fun psi(F: Double, u: Double) = F - m * c * u / 2


}

fun printResult(res : MutableList<MutableList<Double>>) {
    println("╔═════════════╦═════════════╦═════════════╦═════════════╗")
    println("║      Z      ║     F(z)    ║     u(z)    ║    u_p(z)   ║")
    println("╠═════════════╬═════════════╬═════════════╬═════════════╣")
    for (i in 0 until res[0].size) {
        println("║  %8.3f   ║  %8.3f   ║  %8.3e  ║  %8.3e  ║".format(
            res[0][i], res[1][i], res[2][i], res[3][i]))
        if (i != res[0].size - 1) {
            println("╠═════════════╬═════════════╬═════════════╬═════════════╣")
        }
        else
            println("╚═════════════╩═════════════╩═════════════╩═════════════╝")
    }
}

fun main(args: Array<String>) {
    val solver = RungeKuttSolver(0.35, 0.786, 0.0008, 10000.0, 2000.0, 4)
    val n = 30.0
    val z_start = 0.0
    val z_stop = 1.0
    val f_start = 0.0
    val up = solver.u_p(z_start)
    val eps = 1e-7
    var ksi_l = 0.0
    var ksi_r = 1.0
    while (abs(ksi_l - ksi_r) > eps) {
        val ksi_c = (ksi_l + ksi_r) / 2
        val r_c = solver.runge_4th(z_start, z_stop, n, f_start, ksi_c * up)
        val r_l = solver.runge_4th(z_start, z_stop, n, f_start, ksi_l * up)

        val psi_c = solver.psi(r_c[1][r_c[1].size - 1], r_c[2][r_c[2].size - 1])
        val psi_l = solver.psi(r_l[1][r_l[1].size - 1], r_l[2][r_l[2].size - 1])

        if (psi_c * psi_l < 0)
            ksi_r = ksi_c
        else
            ksi_l = ksi_c
    }
    val ksi = (ksi_l + ksi_r) / 2
//    ksi = 0.15

    println("Оптимальное значение ksi: ${ksi}")
    val res = solver.runge_4th(z_start, z_stop, n, f_start, ksi * up)
    res.add(mutableListOf<Double>())
    for (z in res[0])
        res[3].add(solver.u_p(z))

    printResult(res)
    val plot1 = Plot2DPanel()
    plot1.preferredSize = Dimension(600, 1000)
    plot1.addLinePlot("F(z)", Color.CYAN, res[0].toDoubleArray(), res[1].toDoubleArray())
    plot1.setAxisLabel(0, "Z")
    plot1.setAxisLabel(1, "F(Z), Вт/см^2")
    val plot2 = Plot2DPanel()
    plot2.preferredSize = Dimension(600, 1000)
    plot2.addLinePlot("u(z)", Color.CYAN, res[0].toDoubleArray(), res[2].toDoubleArray())
    plot2.setAxisLabel(0, "Z")
    plot2.setAxisLabel(1, "u(Z), Дж/см^3")
    val plot3 = Plot2DPanel()
    plot3.preferredSize = Dimension(600, 1000)
    plot3.addLinePlot("u_p(z)", Color.CYAN, res[0].toDoubleArray(), res[3].toDoubleArray())
    plot3.setAxisLabel(0, "Z")
    plot3.setAxisLabel(1, "u_p(Z)")
    val frame = JFrame()
    frame.extendedState = JFrame.MAXIMIZED_BOTH
//    frame.layout = GridLayout()
//    frame.add(panel1, 0)
//    frame.add(panel2, 1)
    frame.add(plot1, BorderLayout.EAST)
    frame.add(plot2, BorderLayout.CENTER)
    frame.add(plot3, BorderLayout.WEST)


    frame.isVisible = true
    frame.defaultCloseOperation = 2

    println("F0 = ${res[1][0]}\n f1 = ${res[1][res[1].size - 1]}\nu0= ${res[2][0]}\nu1 = ${res[2][res[2].size - 1]}")

}
