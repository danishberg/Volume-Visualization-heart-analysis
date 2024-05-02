// Tutorial: Command line: Step 1 - javac CW.java, Step 2 - java CW 429 180 or java CW.java 429 180. . 
// Important: 429 resolution is considered maximum for my pc, if higher you could encounter this:
// If Error: Exception in thread "main" java.lang.OutOfMemoryError: Java heap space at CW.sampleHeartEquation(CW.java:27) 
// Recommended range: from min: 300 to max: 429      

import javax.imageio.ImageIO;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

public class CW {

    public static void main(String[] args) {
        // Command line parsing and input validation
        if (args.length != 2) {
            System.err.println("Usage: java CW <grid resolution> <rotation angle>");
            return;
        }
        int resolution = Integer.parseInt(args[0]);
        double rotationAngleDegrees = Double.parseDouble(args[1]); // Parse the rotation angle
    
        double[][][] grid = sampleHeartEquation(resolution);
        Vector3[][][] gradients = calculateGradients(grid, resolution);
        BufferedImage image = raycastAndShade(grid, gradients, resolution, resolution, resolution, rotationAngleDegrees);
        saveImage(image, "result.tiff");
    }

     // Sampling the heart equation within a 3D grid
    public static double[][][] sampleHeartEquation(int resolution) {
        double[][][] grid = new double[resolution][resolution][resolution];
        double step = 4.0 / (resolution - 1);
        for (int i = 0; i < resolution; i++) {
            for (int j = 0; j < resolution; j++) {
                for (int k = 0; k < resolution; k++) {
                    double x = i * step - 2;
                    double y = j * step - 2;
                    double z = k * step - 2;
                    grid[i][j][k] = heartEquation(x, y, z);
                }
            }
        }
        return grid;
    }

    // Heart equation definition // Self-Note: DO NOT CHANGE AT ALL 
    public static double heartEquation(double x, double y, double z) {
        return -(Math.pow(x * x + 2 * y * y + z * z - 1, 3) - x * x * z * z * z - 0.1 * y * y * z * z * z);
    }

    // Calculating gradient vectors at each sampling location on the grid
    public static Vector3[][][] calculateGradients(double[][][] grid, int resolution) {
        Vector3[][][] gradients = new Vector3[resolution][resolution][resolution];
        double step = 4.0 / (resolution - 1);
        for (int i = 1; i < resolution - 1; i++) {
            for (int j = 1; j < resolution - 1; j++) {
                for (int k = 1; k < resolution - 1; k++) {
                    double gx = (grid[i + 1][j][k] - grid[i - 1][j][k]) / (2 * step);
                    double gy = (grid[i][j + 1][k] - grid[i][j - 1][k]) / (2 * step);
                    double gz = (grid[i][j][k + 1] - grid[i][j][k - 1]) / (2 * step);
                    gradients[i][j][k] = new Vector3(gx, gy, gz);
                }
            }
        }
        return gradients;
    }

    // Raycasting and shading method with interpolation
    public static BufferedImage raycastAndShade(double[][][] grid, Vector3[][][] gradients, int resolution, int width, int height, double rotationAngleDegrees) {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        double stepX = (double) resolution / width;
        double stepZ = (double) resolution / height;
        double angleRadians = Math.toRadians(rotationAngleDegrees);
        double cosTheta = Math.cos(angleRadians);
        double sinTheta = Math.sin(angleRadians);
    
        for (int px = 0; px < width; px++) {
            for (int pz = 0; pz < height; pz++) {
                boolean hit = false;
                double xWorld = (px - width / 2.0) * stepX;
                double zWorld = (pz - height / 2.0) * stepZ;
    
                // Rotate around the Y-axis
                double rotatedX = cosTheta * xWorld + sinTheta * zWorld;
                double rotatedZ = -sinTheta * xWorld + cosTheta * zWorld;
    
                int igx = Math.min(Math.max((int)((rotatedX + resolution / 2.0) / stepX), 0), resolution - 2);
                int igz = Math.min(Math.max((int)((rotatedZ + resolution / 2.0) / stepZ), 0), resolution - 2);
    
                for (int gy = 0; gy < resolution - 2; gy++) {
                    if (grid[igx][gy][igz] * grid[igx][gy + 1][igz] <= 0) {
                        // Trilinear interpolation to find the exact surface crossing
                        double t = -grid[igx][gy][igz] / (grid[igx][gy + 1][igz] - grid[igx][gy][igz]);
                        Vector3 grad1 = gradients[igx][gy][igz];
                        Vector3 grad2 = gradients[igx][gy + 1][igz];
                        Vector3 normal = interpolateNormal(grad1, grad2, t).normalize();
    
                        // Calculate lighting // Note: Changing Ambient shows better color intensity
                        double ambient = 0.5;
                        double diffuse = 0.7 * Math.max(0, Vector3.dotProduct(normal, new Vector3(1, -1, 1).normalize()));
                        double specular = 0.3 * Math.pow(Math.max(0, Vector3.dotProduct(new Vector3(0, 0, 1), normal.reflect(new Vector3(1, -1, 1).normalize()))), 20);
    
                        double brightness = ambient + diffuse + specular;
                        int colorValue = (int) (brightness * 255);
                        colorValue = Math.min(colorValue, 255);
                        int color = (colorValue << 16) | (colorValue << 8) | colorValue;
                        image.setRGB(px, pz, color);
                        hit = true;
                        break;
                    }
                }
    
                if (!hit) {
                    image.setRGB(px, pz, Color.BLACK.getRGB());
                }
            }
        }
    
        return image;
    }
    
    // Trilinear Interpolation helper method for normals between two gradients based on t
    private static Vector3 interpolateNormal(Vector3 grad1, Vector3 grad2, double t) {
        double x = grad1.x * (1 - t) + grad2.x * t;
        double y = grad1.y * (1 - t) + grad2.y * t;
        double z = grad1.z * (1 - t) + grad2.z * t;
        return new Vector3(x, y, z);
    }
   // Saving the final image to reqired TIFF format file, PNG could be used if both 21 and 125 code lines are changed accordingly
    public static void saveImage(BufferedImage image, String filename) {
        try {
            File outputFile = new File(filename);
            ImageIO.write(image, "tiff", outputFile);
            System.out.println("Image saved as " + filename);
        } catch (IOException e) {
            System.err.println("Error saving image: " + e.getMessage());
            e.printStackTrace();
        }
    }
}

// Vector3 helper class for gradient calculations and basic vector math
class Vector3 {
    public double x, y, z;

    public Vector3(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public Vector3 add(Vector3 other) {
        return new Vector3(x + other.x, y + other.y, z + other.z);
    }

    public Vector3 subtract(Vector3 other) {
        return new Vector3(x - other.x, y - other.y, z - other.z);
    }

    public Vector3 multiply(double scalar) {
        return new Vector3(x * scalar, y * scalar, z * scalar);
    }

    public double length() {
        return Math.sqrt(x * x + y * y + z * z);
    }

    public Vector3 normalize() {
        double len = length();
        if (len == 0) return new Vector3(0, 0, 0);
        return new Vector3(x / len, y / len, z / len);
    }
    
    public Vector3 reflect(Vector3 normal) {
        double dotProduct = Vector3.dotProduct(this, normal);
        return normal.multiply(2 * dotProduct).subtract(this);
    }

    public static double dotProduct(Vector3 a, Vector3 b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
}
