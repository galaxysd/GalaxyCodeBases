/*
 DNA Lamp
 Hardwre: https://store.arduino.cc/usa/arduino-nano
 Inspired: http://www.instructables.com/id/3D-Printed-DNA-Lamp/
 Morse Code: http://www.instructables.com/id/Arduino-Morse-Code/
 Aonther better Morse Code: https://create.arduino.cc/projecthub/team-onyx/morse-code-translator-46e9b8
*/

int ledB = 10;         // the pin that the Base LEDs are attached to
int led = 9;           // the pin that the Top LEDs are attached to
int button = 8;        // the pin that the button is connected to
int motor = 11;        // the pin that the relay for the motor is attached to
volatile int mode=0;            // initial operation mode

int twinkle=36;
int ttime=3;
// http://www.nu-ware.com/NuCode%20Help/index.html?morse_code_structure_and_timing_.htm
int audioO = 6;      // output audio on pin 8
int note = 1200;      // music note/pitch
// const char stringToMorseCode[] = "SOS  .  This  is  Galaxy  @  BGI  -  SZ  .  ";
const char stringToMorseCode[] = "SOS  ";
int dotLen = 150;     // length of the morse code 'dot'
int dashLen = dotLen * 3;    // length of the morse code 'dash'
int elemPause = dotLen;  // length of the pause between elements of a character
int Spaces = dotLen * 3;     // length of the spaces between characters
int wordPause = dotLen * 7;  // length of the pause between words
int charPause = dotLen * 2;  // length of the pause between characters (3-1=2)

int MorseHigh = 255;
int MorseLow = 8;
volatile int MorseLight = MorseLow;

int brightness = 0;    // how bright the LED is
int fadeAmount = 5;    // how many points to fade the LED by

// the setup routine runs once when you press reset:
void setup()  {
  pinMode(led, OUTPUT);          // declare pin 9 to be an output:
  pinMode(ledB, OUTPUT);          // declare pin 9 to be an output:
  pinMode(motor, OUTPUT);        // declare pin 11 to be an output:
  pinMode(button, INPUT_PULLUP); // declare pin 8 to be an input with internal pullups:
}

// the loop routine runs over and over again forever:
void loop()  {

  getBtn();

  if (mode==0) //LEDs fading and rotation motor ON
  {
    // set the brightness of pin 9:
    analogWrite(led, brightness);
    analogWrite(ledB, brightness^0xFF);

    if (brightness==0 or brightness==255)
    {
      delay(500); //time while fading that LEDs will be off
    }

    // change the brightness for next time through the loop:
    brightness = brightness + fadeAmount;

    // reverse the direction of the fading at the ends of the fade:
    if (brightness == 0 || brightness == 255) {
      fadeAmount = -fadeAmount ;
    }
    // wait for 30 milliseconds to see the dimming effect
    delay(30);
    digitalWrite(motor, 1);
  }

  else if (mode==1) //LEDs fading and rotation motor OFF
  {
    // set the brightness of pin 9:
    analogWrite(led, brightness);
    analogWrite(ledB, brightness^0xFF);

    if (brightness==0 or brightness==255)
    {
      delay(1000);
    }

    // change the brightness for next time through the loop:
    brightness = brightness + fadeAmount;

    // reverse the direction of the fading at the ends of the fade:
    if (brightness == 0 || brightness == 255) {
      fadeAmount = -fadeAmount ;
    }
    // wait for 30 milliseconds to see the dimming effect
    delay(30);
    digitalWrite(motor, 0);
  }

  else if (mode==2) //LEDs MorseLow and rotation motor ON
  {
    digitalWrite(motor, 1);
    digitalWrite(led, 0);
    digitalWrite(ledB, 0);
    delay(50);
	MorseLight = MorseLow;
    for (int i = 0; i < sizeof(stringToMorseCode) - 1; i++) {
        // Get the character in the current position
        char tmpChar = stringToMorseCode[i];
        // Set the case to lower case
        tmpChar = toLowerCase(tmpChar);
        // Call the subroutine to get the morse code equivalent for this character
        GetChar(tmpChar);
        delay(charPause);
    }
  }

  else if (mode==3) //LEDs ON and rotation motor OFF
  {
    digitalWrite(motor, 0);
    digitalWrite(led, 1);
    digitalWrite(ledB, 1);
    delay(50);
  }

  else if (mode==4) //LEDs MorseHigh and rotation motor OFF
  {
    digitalWrite(motor, 0);
    digitalWrite(led, 0);
    digitalWrite(ledB, 0);
    delay(50);
	MorseLight = MorseHigh;
    for (int i = 0; i < sizeof(stringToMorseCode) - 1; i++) {
        // Get the character in the current position
        char tmpChar = stringToMorseCode[i];
        // Set the case to lower case
        tmpChar = toLowerCase(tmpChar);
        // Call the subroutine to get the morse code equivalent for this character
        GetChar(tmpChar);
		delay(charPause);
    }
	//delay(wordPause);
  }

  else if (mode==5) //LEDs OFF and rotation motor OFF
  {
    digitalWrite(motor, 0);
    digitalWrite(led, 0);
    digitalWrite(ledB, 0);
    delay(50);
  }
}

// DOT
void MorseDot()
{
  analogWrite(ledB, MorseLight);  	// turn the LED on 
  analogWrite(led, MorseLight);
  tone(audioO, note, dotLen);	// start playing a tone
  delay(dotLen);             	// hold in this position
}

// DASH
void MorseDash()
{
  analogWrite(ledB, MorseLight);  	// turn the LED on 
  analogWrite(led, MorseLight);
  tone(audioO, note, dashLen);	// start playing a tone
  delay(dashLen);               // hold in this position
}

// Turn Off
void LightsOff(int delayTime)
{
  digitalWrite(ledB, LOW);    	// turn the LED off  	
  digitalWrite(led, LOW);
  noTone(audioO);	       	   	// stop playing a tone
  delay(delayTime);            	// hold in this position
}

void getBtn() {
    if (digitalRead(button)==LOW) //changes operation mode when switch is pressed
    {
      if (mode!=5) {
        mode=mode+1;
      } else {
        mode=0;
      }
      for(int x=0; x<ttime; x++) {
          analogWrite(ledB, MorseLow); delay(twinkle);digitalWrite(ledB, 0); delay(twinkle);
      }
      analogWrite(ledB, brightness^0xFF);
      //delay(abs(1500-2*ttime*twinkle)); //time to wait between push button can be pressed again
    }
}

// *** Characters to Morse Code Conversion *** //
void GetChar(char tmpChar) {
	getBtn();
	// Take the passed character and use a switch case to find the morse code for that character. https://morsecode.scphillips.com/morse2.html
	switch (tmpChar) {
	  case 'a':	
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'b':
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'c':
	    MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'd':
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'e':
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'f':
	    MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'g':
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'h':
	    MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'i':
	    MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'j':
	    MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
      case 'k':
	    MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'l':
	    MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
      case 'm':
	    MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'n':
	    MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'o':
	    MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'p':
	    MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 'q':
	    MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'r':
	    MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 's':
	    MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case 't':
	    MorseDash();
		LightsOff(elemPause);
		break;
	  case 'u':
	    MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'v':
	    MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'w':
	    MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'x':
	    MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'y':
	    MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  case 'z':
	    MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
  	  case '.':
  	    MorseDot();
  		LightsOff(elemPause);
  		MorseDash();
  		LightsOff(elemPause);
  	    MorseDot();
  		LightsOff(elemPause);
  		MorseDash();
  		LightsOff(elemPause);
  	    MorseDot();
  		LightsOff(elemPause);
  		MorseDash();
  		LightsOff(elemPause);
  		break;
  	  case ',':
  		MorseDash();
  		LightsOff(elemPause);
  		MorseDash();
  		LightsOff(elemPause);
  		MorseDot();
  		LightsOff(elemPause);
  		MorseDot();
  		LightsOff(elemPause);
  		MorseDash();
  		LightsOff(elemPause);
  		MorseDash();
  		LightsOff(elemPause);
  		break;
	  case '?':
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case '@': // .--.-.
	    MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
	    MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
	    MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
		break;
	  case '-': // -....-
	    MorseDash();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
	    MorseDot();
		LightsOff(elemPause);
		MorseDot();
		LightsOff(elemPause);
	    MorseDot();
		LightsOff(elemPause);
		MorseDash();
		LightsOff(elemPause);
		break;
	  default: 
		// If a matching character was not found it will default to a blank space
		LightsOff(Spaces);			
	}
}

