import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.*;
import javax.imageio.ImageIO;
import java.util.Scanner;

/**
 *
 * @author Dylan Howey <dhowey@cord.edu>
 * @since  2012-5-2
 */
class PrintArray {

    public static void main(String[] args) {
    	int width = 0;
	int height = 0;
	String filename = "";
	int[][] ferro = new int[1000][1000];

	try {
            width = Integer.parseInt(args[0]); //The number of columns of the input
            height = Integer.parseInt(args[1]); //The number of rows of the input
            filename = args[2]; //The name of the input file
	} catch (Exception e) {
            System.out.println("Usage: Ising.jar [width] [height] [input file]");
            System.exit(0);
        }

	if(width > 1000 || height > 1000) {
		System.out.println("Error: width and height must both be less than 1000");
		System.exit(1);
	}
	
	BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        WritableRaster raster = image.getRaster();

	Scanner in = new Scanner(System.in); //Stupid hack because I'm dumb with Java
	try {
     		in = new Scanner(new FileInputStream(filename));
	} catch (Exception e) {
		System.out.println("Error opening " + filename + " for reading");
		System.exit(1);
	}

	for(int row = 0; row < height; row++) {
		for(int col = 0; col < width; col++) {
			ferro[row][col] = in.nextInt();
		}
	}

        int[] rgb = new int[3]; //Index 0 - red; Index 1 - green; Index 2 - blue
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                //If the dipole is pointing one way, paint it a certain color
                if (ferro[row][col] > 0) {
                    rgb[0] = 255;
                    rgb[1] = 255;
                    rgb[2] = 0;
                } //Else, paint it a different color
                else {
                    rgb[0] = 0;
                    rgb[1] = 0;
                    rgb[2] = 255;

                }
                raster.setPixel(row, col, rgb);
            }
        }
        try {
            //Writes the image to a file
            ImageIO.write(image, "jpg", new File("output.jpg"));
        } catch (IOException e) {
            System.out.println("Error saving file 'output.jpg'");
        }
    }
}
