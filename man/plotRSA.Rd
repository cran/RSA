\name{plotRSA}
\alias{plotRSA}
\title{Plots a response surface of a polynomial equation of second degree}
\usage{
  plotRSA(x = NULL, y = 0, x2 = 0, y2 = 0, xy = 0, w = 0,
    wx = 0, wy = 0, b0 = 0, xlim = NULL, ylim = NULL,
    zlim = NULL, xlab = NULL, ylab = NULL, zlab = NULL,
    surface = "predict", lambda = NULL,
    rot = list(x = -45, y = 45, z = 35),
    label.rot = list(x = 45, y = -25, z = 94),
    gridsize = 21, bw = FALSE, legend = TRUE,
    showSP = TRUE, axes = c("LOC", "LOIC", "PA1", "PA2"),
    cex = 1.2, type = "3d", points = FALSE, model = "full",
    demo = FALSE, fit = NULL, link = "identity",
    distance = c(1, 1, 1), tck = c(1, 1, 1),
    showBorder = TRUE, showContour = FALSE, pal = NULL,
    ...)
}
\arguments{
  \item{x}{Either an RSA object (returned by the \code{RSA}
  function), or the coefficient for the X predictor}

  \item{y}{Y coefficient}

  \item{x2}{X^2 coefficient}

  \item{y2}{Y^2 coefficient}

  \item{xy}{XY interaction coefficient}

  \item{w}{W coefficient (for (un)constrained absolute
  difference model)}

  \item{wx}{WX coefficient (for (un)constrained absolute
  difference model)}

  \item{wy}{WY coefficient (for (un)constrained absolute
  difference model)}

  \item{b0}{Intercept}

  \item{xlim}{Limits of the x axis}

  \item{ylim}{Limits of the y axis}

  \item{zlim}{Limits of the z axis}

  \item{xlab}{Label for x axis}

  \item{ylab}{Label for y axis}

  \item{zlab}{Label for z axis}

  \item{surface}{Method for the calculation of the surface
  z values. "predict" takes the predicted values from the
  model, "smooth" uses a thin plate smoother (function
  \code{Tps} from the \code{fields} package) of the raw
  data}

  \item{lambda}{lambda parameter for the smoother. Default
  (NULL) means that it is estimated by the smoother
  function. Small lambdas around 1 lead to rugged surfaces,
  big lambdas to very smooth surfaces.}

  \item{rot}{Rotation of the 3d surface plot (when type ==
  "3d")}

  \item{label.rot}{Rotation of the axis labls (when type ==
  "3d")}

  \item{gridsize}{Number of grid nodes in each dimension}

  \item{bw}{Print surface in black and white instead of
  colors?}

  \item{legend}{Print color legend for z values?}

  \item{cex}{Font size factor for axes labels}

  \item{type}{\code{3d} for 3d surface plot, \code{contour}
  for 2d contour plot, "interactive" for interactive
  rotatable plot. Shortcuts (i.e., first letter of string)
  are sufficient}

  \item{points}{Should the original data points be
  overplotted?}

  \item{model}{If x is an RSA object: from which model
  should the response surface be computed?}

  \item{demo}{Do not change that parameter (internal use
  only)}

  \item{fit}{Do not change that parameter (internal use
  only)}

  \item{showSP}{Should the surface parameters a1 to a4 be
  shown on the plot? In case of a 3d plot a1 to a4 are
  printed on the upper left side; in case of a contour plot
  the principal axes are plotted.}

  \item{axes}{A vector of strings specifying the axes that
  should be plotted. Can be any combination of c("LOC",
  "LOIC", "PA1", "PA2"). LOC = line of congruence, LOIC =
  line of incongruence, PA1 = first principal axis, PA2 =
  second principal axis}

  \item{link}{Link function to transform the z axes.
  Implemented are "identity" (no transformation; default),
  "probit", and "logit"}

  \item{showBorder}{Should a thicker border around the
  surface be plotted? Sometimes this border leaves the
  surrounding box, which does not look good. In this case
  the border can be suppressed by setting
  \code{showBorder=FALSE}.}

  \item{showContour}{Should the contour lines be plotted on
  the 3d wireframe plot? (Parameter only relevant for
  type="3d")}

  \item{distance}{A vector of three values defining the
  distance of labels to the axes}

  \item{tck}{A vector of three values defining the position
  of labels to the axes (see ?wireframe)}

  \item{pal}{A palette for shading}

  \item{...}{Additional parameters apssed to the plotting
  function (e.g., main="Title"). A useful title might be
  the R squared of the plotted model: main =
  as.expression(bquote(R^2==.(round(getPar(x, "r2",
  model="full"), 3))))}
}
\description{
  Plots an RSA object, or a response surface with specified
  parameters
}
\details{
  Each plot type has its distinctive advantages. The
  two-dimensional contour plot gives a clear view of the
  position of the principal axes and the stationary point.
  The 3d plot gives a three dimensional impression of the
  surface, allows overplotting of the original data points
  (in case an RSA object is provided), and allows the
  interactive adjustment of regression weights in the
  \code{\link{demoRSA}} function. The interactive plot
  allows rotating and exploring a three-dimensional surface
  with the mouse (nice for demonstration purposes). If you
  want to export publication-ready plots, it is recommended
  to export it with following commands: \code{p1 <-
  plot(r1, bw=TRUE) trellis.device(device="cairo_pdf",
  filename="RSA_plot.pdf") print(p1) dev.off()}
}
\examples{
# Plot response surfaces from known parameters
# example of Edwards (2002), Figure 3
plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628)
plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, type="c")

# Plot response surface from an RSA object
set.seed(0xBEEF)
n <- 300
err <- 2
x <- rnorm(n, 0, 5)
y <- rnorm(n, 0, 5)
df <- data.frame(x, y)
df <- within(df, {
	diff <- x-y
	absdiff <- abs(x-y)
	sqdiff <- (x-y)^2
	z.diff <- diff + rnorm(n, 0, err)
	z.abs <- absdiff + rnorm(n, 0, err)
	z.sq <- sqdiff + rnorm(n, 0, err)
	z.add <- diff + 0.4*x + rnorm(n, 0, err)
	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
})

r1 <- RSA(z.sq~x*y, df)
plot(r1)
plot(r1, points=TRUE, model="sqdiff")
}
\seealso{
  \code{\link{demoRSA}}, \code{\link{RSA}}
}

